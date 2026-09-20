"""PySCF screening engine: PBE0 + DKH2 + def2-TZVP, opt + harmonic G(298.15 K).

This is not ORCA PBE0-D4 ZORA. Hg uses the def2 ECP that belongs with def2-TZVP
in PySCF (20 valence electrons). DKH2 is the Wolf/Reiher/Hess scalar 1e
Hamiltonian (even1+even2) built from T, V, pVp; ECP is then added to hcore.
def2-SVP is used only if TZVP runs out of memory.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from .geometry import read_xyz, write_xyz
from .species import SPECIES_META

HARTREE_TO_KCAL = 627.509474
T_K = 298.15
P_PA = 101325
TZVP = "def2-TZVP"
SVP = "def2-SVP"


def method_label(basis: str) -> str:
    return (
        f"PySCF PBE0 DKH2 {basis} (Hg def2-ECP with {basis}; "
        "scalar DKH2 1e + ECP; not ORCA PBE0-D4 ZORA)"
    )


def _is_oom(exc: BaseException) -> bool:
    msg = f"{type(exc).__name__}: {exc}".lower()
    return isinstance(exc, MemoryError) or any(
        k in msg for k in ("memory error", "out of memory", "malloc", "std::bad_alloc", "enomem")
    )


def xyz_to_atom(xyz_path: Path) -> str:
    return "; ".join(f"{el} {x:.10f} {y:.10f} {z:.10f}" for el, x, y, z in read_xyz(xyz_path))


def write_mol_xyz(mol, path: Path, comment: str = "pyscf optimized") -> None:
    coords = mol.atom_coords(unit="Angstrom")
    atoms = [(mol.atom_symbol(i), float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2])) for i in range(mol.natm)]
    write_xyz(path, atoms, comment)


def make_mol(xyz_path: Path, charge: int, mult: int, basis: str, max_memory_mb: int):
    from pyscf import gto

    symbols = {el for el, *_ in read_xyz(xyz_path)}
    ecp = {el: basis for el in symbols if el == "Hg"}
    mol = gto.M(
        atom=xyz_to_atom(xyz_path),
        basis=basis,
        ecp=ecp or None,
        charge=int(charge),
        spin=int(mult) - 1,
        unit="Angstrom",
        verbose=3,
        max_memory=max_memory_mb,
    )
    return mol


def dkh2_hcore(mol) -> np.ndarray:
    """Scalar DKH2 1e Hamiltonian (Wolf, Reiher, Hess, J. Chem. Phys. 117, 9215 (2002))."""
    from pyscf.lib.parameters import LIGHT_SPEED as C

    s = np.asarray(mol.intor_symmetric("int1e_ovlp"), dtype=float)
    t = np.asarray(mol.intor_symmetric("int1e_kin"), dtype=float)
    v = np.asarray(mol.intor_symmetric("int1e_nuc"), dtype=float)
    pvp = np.asarray(mol.intor_symmetric("int1e_pnucp"), dtype=float)
    c2 = float(C) * float(C)

    se, su = np.linalg.eigh(s)
    keep = se > 1e-8
    x = su[:, keep] * (1.0 / np.sqrt(se[keep]))
    torth = x.T @ t @ x
    tt, eig = np.linalg.eigh(torth)
    tt = np.clip(tt, 0.0, None)
    # Stable relativistic kinetic: c^2 (sqrt(1+2T/c^2)-1)
    ev0t = 2.0 * tt / (np.sqrt(1.0 + 2.0 * tt / c2) + 1.0)
    e = ev0t + c2
    aa = np.sqrt((c2 + e) / (2.0 * e))
    rr = float(C) / (c2 + e)
    u = x @ eig
    revt = s @ u
    vt = u.T @ v @ u
    pvpt = u.T @ pvp @ u
    ev1t = (aa[:, None] * vt * aa[None, :]) + (
        aa[:, None] * rr[:, None] * pvpt * rr[None, :] * aa[None, :]
    )
    ev2t = _even2c(vt, pvpt, aa, rr, tt, e)
    hkin = (revt * ev0t) @ revt.T
    h = hkin + revt @ (ev1t + ev2t) @ revt.T
    h = 0.5 * (h + h.T)
    return h


def _mat_axa(p, a):
    return p * a[:, None] * a[None, :]


def _mat_arxra(p, a, r):
    return p * a[:, None] * a[None, :] * r[:, None] * r[None, :]


def _mat_1_over_h(p, e):
    return p / (e[:, None] + e[None, :])


def _mat_mulm(p, q, r, t, rr, alpha, beta):
    qtemp = q * (2.0 * t * rr * rr)[None, :]
    return alpha * (qtemp @ r) + beta * p


def _mat_muld(p, q, r, t, rr, alpha, beta):
    t_safe = np.clip(t, 1e-16, None)
    qtemp = q * (0.5 / (t_safe * rr * rr))[None, :]
    return alpha * (qtemp @ r) + beta * p


def _even2c(vv, gg, aa, rr, tt, e):
    v = _mat_axa(np.array(vv, copy=True), aa)
    pvp = _mat_arxra(np.array(gg, copy=True), aa, rr)
    vh = _mat_axa(_mat_1_over_h(np.array(vv, copy=True), e), aa)
    pvph = _mat_arxra(_mat_1_over_h(np.array(gg, copy=True), e), aa, rr)
    w1o1 = -(pvph @ v)
    w1o1 = _mat_muld(w1o1, pvph, pvp, tt, rr, 1.0, 1.0)
    w1o1 = _mat_mulm(w1o1, vh, v, tt, rr, 1.0, 1.0)
    w1o1 = w1o1 - (vh @ pvp)
    o1w1 = pvp @ vh
    o1w1 = _mat_muld(o1w1, pvp, pvph, tt, rr, -1.0, 1.0)
    o1w1 = _mat_mulm(o1w1, v, vh, tt, rr, -1.0, 1.0)
    o1w1 = o1w1 + (v @ pvph)
    return 0.5 * w1o1 - 0.5 * o1w1


def make_mf(mol, density_fit: bool = True):
    from pyscf import dft

    mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
    mf.xc = "PBE0"
    mf.conv_tol = 1e-8
    mf.max_cycle = 128
    if density_fit:
        mf = mf.density_fit()

    def _hcore(mol_h=None):
        m = mf.mol if mol_h is None else mol_h
        h = dkh2_hcore(m)
        if len(getattr(m, "_ecpbas", [])) > 0:
            h = h + m.intor_symmetric("ECPscalar")
        return h

    mf.get_hcore = _hcore
    return mf


def _blank(name: str, charge: int, mult: int, basis: str) -> dict:
    return {
        "species": name,
        "charge": charge,
        "mult": mult,
        "method": method_label(basis),
        "basis": basis,
        "xc": "PBE0",
        "relativistic": "DKH2",
        "temperature_K": T_K,
        "converged": False,
        "opt_converged": None,
        "freq_ok": False,
        "n_imaginary": None,
        "E_elec_Eh": None,
        "E_opt_Eh": None,
        "G_opt_Eh": None,
        "Gcorr_Eh": None,
        "E_sp_Eh": None,
        "G_composite_Eh": None,
        "G_composite_kcal": None,
        "g_source": None,
        "missing_sp": False,
        "error": None,
        "seconds": None,
        "oom": False,
    }


def _count_imag(freq_au) -> int:
    f = np.asarray(freq_au)
    if f.size == 0:
        return 0
    if np.iscomplexobj(f):
        return int(np.sum(np.abs(f.imag) > 1e-12))
    return int(np.sum(np.real(f) < -1e-8))


def _real_freq_au(freq_au) -> np.ndarray:
    f = np.asarray(freq_au, dtype=complex)
    out = np.real(f).astype(float)
    out = np.where(np.abs(f.imag) > 1e-12, 0.0, out)
    out = np.where(out < 0.0, 0.0, out)
    return out


def _ensure_scratch(jobdir: Path) -> Path:
    # WSL /tmp is often a small tmpfs; DF Hessian needs several GB on disk.
    scratch = Path.home() / ".cache" / "pyscf_scratch" / jobdir.name
    scratch.mkdir(parents=True, exist_ok=True)
    s = str(scratch)
    os.environ["TMPDIR"] = s
    os.environ["TMP"] = s
    os.environ["TEMP"] = s
    os.environ["PYSCF_TMPDIR"] = s
    try:
        from pyscf import lib
        lib.param.TMPDIR = s
    except Exception:
        pass
    return scratch


def _pbe0_mf(mol, density_fit: bool = True):
    from pyscf import dft

    mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
    mf.xc = "PBE0"
    mf.conv_tol = 1e-8
    mf.max_cycle = 128
    if density_fit:
        mf = mf.density_fit()
    return mf


def _fd_hessian_from_grad(mf, delta=0.005):
    """Central-difference Hessian from analytic gradients. delta in Bohr."""
    import gc

    mol = mf.mol
    natm = mol.natm
    coords0 = np.asarray(mol.atom_coords(), dtype=float)
    h = np.zeros((natm, 3, natm, 3))
    use_df = bool(getattr(mf, "with_df", None))

    def grad_at(coords):
        mol1 = mol.set_geom_(coords, unit="Bohr", inplace=False)
        mf1 = _pbe0_mf(mol1, density_fit=use_df)
        e = float(mf1.kernel())
        if not bool(getattr(mf1, "converged", False)):
            raise RuntimeError(f"FD Hessian SCF not converged E={e}")
        g = np.asarray(mf1.nuc_grad_method().kernel(), dtype=float)
        return g

    for i in range(natm):
        for x in range(3):
            cp = coords0.copy()
            cm = coords0.copy()
            cp[i, x] += delta
            cm[i, x] -= delta
            gp = grad_at(cp)
            gm = grad_at(cm)
            h[i, x] = (gp - gm) / (2.0 * delta)
            print(f"FD Hessian atom {i} coord {x} done", flush=True)
            gc.collect()
    h = 0.5 * (h + h.transpose(2, 3, 0, 1))
    # PySCF thermo.harmonic_analysis expects (natm, natm, 3, 3)
    return np.transpose(h, (0, 2, 1, 3))


def _harmonic_G(mf, log_path: Path) -> tuple[float | None, float | None, int | None, str | None, str]:
    """G = E_DKH2 + (G_PBE0_harmonic - E_PBE0). Return (G, Gcorr, n_imag, error, g_source)."""
    from pyscf.geomopt.geometric_solver import optimize
    from pyscf.hessian import thermo as th

    e_dkh2 = float(mf.e_tot)
    natm = mf.mol.natm
    if natm == 1:
        info = th.thermo(mf, np.zeros(0), T_K, P_PA)
        g = float(info["G_tot"][0])
        return g, g - e_dkh2, 0, None, "dkh2_elec+atomic_translational_G"
    last_err = None
    mf_h = None
    hessian = None
    e_pbe0 = None
    used_fd = False
    hess_note = "PBE0/def2-TZVP DF harmonic at PBE0 stationary point (not DKH2 picture-change Hessian)"
    # Full (non-DF) hybrid Hessian OOMs on ~8 GB; stay on DF and reopt on the
    # PBE0 surface so frequencies are not evaluated at a non-stationary point.
    try:
        print("PBE0 DF reopt before Hessian", flush=True)
        mf_opt = _pbe0_mf(mf.mol, density_fit=True)
        mol_h = optimize(mf_opt, maxsteps=25)
        write_mol_xyz(mol_h, log_path.parent / "opt_pbe0_hess.xyz", "PBE0/def2-TZVP Hessian geometry")
        mf_h = _pbe0_mf(mol_h, density_fit=True)
        e_pbe0 = float(mf_h.kernel())
        if not bool(getattr(mf_h, "converged", False)):
            return None, None, None, "PBE0 Hessian SCF not converged", "failed"
        n_vib = max(1, 3 * natm - 6)
        # Analytic DF hybrid Hessians are often all-imaginary for anions (and the
        # full non-DF Hessian OOMs). Use FD gradient Hessian when |charge|>=2.
        if abs(int(mf.mol.charge)) >= 2:
            print("FD gradient Hessian (|charge|>=2; skip analytic DF Hessian)", flush=True)
            hessian = _fd_hessian_from_grad(mf_h)
            used_fd = True
            hess_note = "PBE0/def2-TZVP DF finite-difference gradient Hessian (not DKH2 picture-change)"
        else:
            print("PBE0 DF analytic Hessian", flush=True)
            hessian = mf_h.Hessian().kernel()
            freq_info = th.harmonic_analysis(mf_h.mol, hessian, imaginary_freq=True)
            n_imag = _count_imag(freq_info["freq_au"])
            if n_imag >= n_vib:
                last_err = f"DF analytic Hessian: {n_imag}/{n_vib} imaginary; FD gradient Hessian"
                print(last_err, flush=True)
                hessian = _fd_hessian_from_grad(mf_h)
                used_fd = True
                hess_note = "PBE0/def2-TZVP DF finite-difference gradient Hessian (not DKH2 picture-change)"
    except Exception as exc:
        last_err = f"Hessian failed: {type(exc).__name__}: {exc}"
        print(last_err, flush=True)
        if used_fd or mf_h is None or not bool(getattr(mf_h, "converged", False)):
            return None, None, None, last_err, "failed"
        try:
            hessian = _fd_hessian_from_grad(mf_h)
            used_fd = True
            hess_note = "PBE0/def2-TZVP DF finite-difference gradient Hessian (not DKH2 picture-change)"
        except Exception as exc2:
            return None, None, None, f"{last_err} | FD fallback: {type(exc2).__name__}: {exc2}", "failed"
    if hessian is None or mf_h is None or e_pbe0 is None:
        return None, None, None, last_err or "Hessian failed", "failed"
    try:
        np.save(log_path.parent / "hessian_pbe0.npy", np.asarray(hessian, dtype=float))
    except Exception:
        pass
    freq_info = th.harmonic_analysis(mf_h.mol, hessian, imaginary_freq=True)
    n_imag = _count_imag(freq_info["freq_au"])
    freq_au = _real_freq_au(freq_info["freq_au"])
    info = th.thermo(mf_h, freq_au, T_K, P_PA)
    g_pbe0 = float(info["G_tot"][0])
    gcorr = g_pbe0 - e_pbe0
    g = e_dkh2 + gcorr
    (log_path.parent / "freq.json").write_text(
        json.dumps(
            {
                "n_imaginary": n_imag,
                "E_pbe0_Eh": e_pbe0,
                "G_pbe0_Eh": g_pbe0,
                "Gcorr_Eh": gcorr,
                "ZPE": info.get("ZPE"),
                "freq_wavenumber": np.asarray(freq_info.get("freq_wavenumber", freq_au)).real.tolist(),
                "hessian_method": hess_note,
            },
            indent=2,
            default=str,
        )
    )
    return g, gcorr, n_imag, None, "E_DKH2 + Gcorr_PBE0_def2TZVP_harmonic"


def run_one(name: str, xyz_path: Path, jobdir: Path, basis: str = TZVP, max_memory_mb: int = 6000, reuse_opt: bool = False) -> dict:
    jobdir = Path(jobdir)
    jobdir.mkdir(parents=True, exist_ok=True)
    charge, mult = SPECIES_META[name]
    rec = _blank(name, charge, mult, basis)
    log_path = jobdir / "run.log"
    t0 = time.time()
    local_xyz = jobdir / f"{name}.xyz"
    shutil.copy2(xyz_path, local_xyz)
    scratch = _ensure_scratch(jobdir)
    old_cwd = os.getcwd()
    try:
        os.chdir(scratch)
    except Exception:
        old_cwd = None
    log_f = log_path.open("w", encoding="utf-8", buffering=1)
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout = sys.stderr = log_f
        n_atom = len(read_xyz(local_xyz))
        opt_xyz = jobdir / "opt.xyz"
        if n_atom == 1:
            mol = make_mol(local_xyz, charge, mult, basis, max_memory_mb)
            rec["opt_converged"] = True
            write_mol_xyz(mol, opt_xyz, f"{name} atom")
        elif reuse_opt and opt_xyz.exists():
            mol = make_mol(opt_xyz, charge, mult, basis, max_memory_mb)
            rec["opt_converged"] = True
        else:
            mol = make_mol(local_xyz, charge, mult, basis, max_memory_mb)
            mf_opt = make_mf(mol, density_fit=True)
            from pyscf.geomopt.geometric_solver import optimize

            try:
                mol_eq = optimize(mf_opt, maxsteps=40)
                rec["opt_converged"] = True
                write_mol_xyz(mol_eq, opt_xyz, f"{name} PBE0 DKH2 {basis} optimized")
                mol = make_mol(opt_xyz, charge, mult, basis, max_memory_mb)
            except Exception as exc:
                rec["opt_converged"] = False
                rec["error"] = f"geometry opt failed: {type(exc).__name__}: {exc}"
                write_mol_xyz(mol, opt_xyz, f"{name} unoptimized starting geometry")
                log_f.write(traceback.format_exc() + "\n")
                mol = make_mol(opt_xyz, charge, mult, basis, max_memory_mb)
        mf = make_mf(mol, density_fit=True)

        e = mf.kernel()
        scf_ok = bool(getattr(mf, "converged", False))
        rec["E_elec_Eh"] = float(e)
        rec["E_sp_Eh"] = float(e)
        rec["E_opt_Eh"] = float(e)
        write_energy_json(jobdir, rec)
        if not scf_ok:
            rec["error"] = "SCF not converged"
            rec["converged"] = False
        else:
            g, gcorr, n_imag, herr, g_source = _harmonic_G(mf, log_path)
            rec["n_imaginary"] = n_imag
            rec["g_source"] = g_source
            if herr:
                rec["error"] = herr
                rec["freq_ok"] = False
                rec["converged"] = False
            else:
                rec["Gcorr_Eh"] = gcorr
                rec["G_opt_Eh"] = g
                rec["G_composite_Eh"] = g
                rec["G_composite_kcal"] = g * HARTREE_TO_KCAL
                rec["missing_sp"] = False
                rec["freq_ok"] = True
                rec["converged"] = bool(rec.get("opt_converged"))
                notes = []
                if rec.get("error"):
                    notes.append(rec["error"])
                if n_imag and n_imag > 0:
                    notes.append(f"{n_imag} imaginary frequency(ies) after harmonic analysis")
                rec["error"] = " | ".join(notes) if notes else None
        (jobdir / "energy.txt").write_text(
            f"E_elec_Eh={rec.get('E_elec_Eh')}\nG_composite_kcal={rec.get('G_composite_kcal')}\n"
            f"converged={rec.get('converged')}\n{rec.get('method')}\n"
        )
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["converged"] = False
        rec["oom"] = _is_oom(exc)
        log_f.write(traceback.format_exc() + "\n")
    finally:
        sys.stdout, sys.stderr = old_out, old_err
        log_f.close()
        if old_cwd:
            try:
                os.chdir(old_cwd)
            except Exception:
                pass
    rec["seconds"] = round(time.time() - t0, 3)
    write_energy_json(jobdir, rec)
    return rec


def write_energy_json(jobdir: Path, rec: dict) -> Path:
    path = Path(jobdir) / "energy.json"
    dump = {k: v for k, v in rec.items() if k != "oom"}
    path.write_text(json.dumps(dump, indent=2, default=str) + "\n")
    return path


def run_species(name: str, xyz_path: Path, jobdir: Path, max_memory_mb: int = 6000, reuse_opt: bool = False) -> dict:
    rec = run_one(name, xyz_path, jobdir, basis=TZVP, max_memory_mb=max_memory_mb, reuse_opt=reuse_opt)
    if rec.get("oom"):
        rec_svp = run_one(name, xyz_path, jobdir, basis=SVP, max_memory_mb=max_memory_mb, reuse_opt=reuse_opt)
        rec_svp["error"] = (
            (rec.get("error") or "TZVP OOM") + " | fallback def2-SVP: " + (rec_svp.get("error") or "ok")
        )
        return rec_svp
    return rec


def has_pyscf() -> bool:
    try:
        import pyscf  # noqa: F401

        return True
    except Exception:
        return False
