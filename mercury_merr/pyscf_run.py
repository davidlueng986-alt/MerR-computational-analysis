from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
import time
import traceback
from pathlib import Path

from .geometry import read_xyz
from .species import SPECIES_META

HARTREE_TO_KCAL = 627.509474
METHOD_LABEL = "PySCF PBE0 def2-TZVP C-PCM(Water), Hg def2-ECP (scalar relativity in ECP; DKH2 not available in PySCF 2.14; not ORCA PBE0-D4 ZORA)"


def win_to_wsl(path: Path) -> str:
    p = str(path.resolve())
    if len(p) >= 2 and p[1] == ":":
        return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")
    return p.replace("\\", "/")


def has_pyscf() -> bool:
    try:
        import pyscf  # noqa: F401
        return True
    except Exception:
        return False


def xyz_to_atom(xyz_path: Path) -> str:
    return "; ".join(f"{el} {x:.10f} {y:.10f} {z:.10f}" for el, x, y, z in read_xyz(xyz_path))


def make_mol(xyz_path: Path, charge: int, mult: int, basis: str, ecp: str, max_memory_mb: int):
    from pyscf import gto
    spin = int(mult) - 1
    mol = gto.M(
        atom=xyz_to_atom(xyz_path),
        basis=basis,
        ecp=ecp,
        charge=int(charge),
        spin=spin,
        unit="Angstrom",
        verbose=3,
        max_memory=max_memory_mb,
    )
    return mol


def make_mf(mol, xc: str, eps: float | None, density_fit: bool):
    from pyscf import dft
    mf = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
    mf.xc = xc
    mf.conv_tol = 1e-8
    mf.max_cycle = 100
    if density_fit:
        mf = mf.density_fit()
    if eps is not None:
        try:
            mf = mf.PCM()
            mf.with_solvent.method = "C-PCM"
            mf.with_solvent.eps = float(eps)
        except Exception:
            pass
    return mf


def run_one(name: str, xyz_path: Path, jobdir: Path, cfg: dict, do_opt: bool = False) -> dict:
    jobdir.mkdir(parents=True, exist_ok=True)
    charge, mult = SPECIES_META[name]
    pc = cfg.get("pyscf", {})
    xc = pc.get("xc", "PBE0")
    basis = pc.get("basis", "def2-tzvp")
    ecp = pc.get("ecp", "def2-tzvp")
    eps = float(pc.get("eps", 78.3553)) if pc.get("solvent", "C-PCM") else None
    density_fit = bool(pc.get("density_fit", True))
    max_memory_mb = int(pc.get("max_memory_mb", 4000))
    do_opt = bool(do_opt or pc.get("optimize", False))
    rec = {
        "species": name,
        "charge": charge,
        "mult": mult,
        "method": METHOD_LABEL,
        "basis": basis,
        "xc": xc,
        "converged": False,
        "opt_converged": None,
        "E_opt_Eh": None,
        "G_opt_Eh": None,
        "Gcorr_Eh": 0.0,
        "E_sp_Eh": None,
        "G_composite_Eh": None,
        "G_composite_kcal": None,
        "error": None,
        "seconds": None,
    }
    log_path = jobdir / "run.log"
    t0 = time.time()
    try:
        mol = make_mol(xyz_path, charge, mult, basis, ecp, max_memory_mb)
        n_atom = mol.natm
        if do_opt and n_atom > 1:
            from pyscf.geomopt.geometric_solver import optimize
            mf_opt = make_mf(mol, xc, eps, density_fit)
            mol_eq = optimize(mf_opt, maxsteps=int(pc.get("opt_maxsteps", 25)))
            e_opt = float(mf_opt.e_tot) if getattr(mf_opt, "e_tot", None) is not None else None
            rec["opt_converged"] = bool(getattr(mf_opt, "converged", False))
            rec["E_opt_Eh"] = e_opt
            xyz_eq = jobdir / "opt.xyz"
            write_mol_xyz(mol_eq, xyz_eq)
            mol = make_mol(xyz_eq, charge, mult, basis, ecp, max_memory_mb)
        mf = make_mf(mol, xc, eps, density_fit)
        e = mf.kernel()
        rec["converged"] = bool(mf.converged)
        rec["E_sp_Eh"] = float(e)
        rec["G_composite_Eh"] = float(e)
        rec["G_composite_kcal"] = float(e) * HARTREE_TO_KCAL
        (jobdir / "energy.txt").write_text(f"{e:.12f} Eh  converged={mf.converged}\n{METHOD_LABEL}\n")
        if not mf.converged:
            rec["error"] = "SCF not converged"
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {exc}"
        rec["converged"] = False
        log_path.write_text(traceback.format_exc())
    rec["seconds"] = round(time.time() - t0, 3)
    (jobdir / "result.json").write_text(json.dumps(rec, indent=2))
    if rec["error"] and not log_path.exists():
        log_path.write_text(rec["error"] + "\n")
    return rec


def _csv_cell(v) -> str:
    if v is None:
        return ""
    s = str(v)
    if any(ch in s for ch in [",", '"', "\n"]):
        return '"' + s.replace('"', '""') + '"'
    return s


def write_mol_xyz(mol, path: Path) -> None:
    coords = mol.atom_coords(unit="Angstrom")
    lines = [str(mol.natm), "pyscf optimized"]
    for i in range(mol.natm):
        x, y, z = coords[i]
        lines.append(f"{mol.atom_symbol(i):2s} {x: .10f} {y: .10f} {z: .10f}")
    path.write_text("\n".join(lines) + "\n")


def write_csv(rows: list[dict], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "species", "E_opt_Eh", "G_opt_Eh", "Gcorr_Eh", "E_sp_Eh",
        "G_composite_Eh", "G_composite_kcal", "converged", "error", "method",
    ]
    lines = [",".join(cols)]
    for r in rows:
        lines.append(",".join(_csv_cell(r.get(c, "")) for c in cols))
    out_csv.write_text("\n".join(lines) + "\n")


def run_all(structures_dir, jobs_dir, out_csv, config_path, do_opt: bool = False) -> list[dict]:
    import yaml
    cfg = yaml.safe_load(Path(config_path).read_text()) if Path(config_path).exists() else {}
    structures_dir = Path(structures_dir)
    jobs_dir = Path(jobs_dir)
    jobs_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in SPECIES_META:
        xyz = structures_dir / f"{name}.xyz"
        rec = {
            "species": name,
            "converged": False,
            "error": f"missing xyz: {xyz}",
            "method": METHOD_LABEL,
            "Gcorr_Eh": 0.0,
        }
        if xyz.exists():
            rec = run_one(name, xyz, jobs_dir / name, cfg, do_opt=do_opt)
        rows.append(rec)
        print(json.dumps({k: rec.get(k) for k in ("species", "converged", "E_sp_Eh", "error", "seconds")}))
    write_csv(rows, Path(out_csv))
    meta = {
        "engine": "pyscf",
        "method": METHOD_LABEL,
        "energies": str(Path(out_csv)),
        "n_converged": sum(1 for r in rows if r.get("converged")),
        "n_species": len(rows),
    }
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    (Path(out_csv).parent / "qm_engine.json").write_text(json.dumps(meta, indent=2))
    return rows


def reexec_wsl(argv: list[str]) -> int:
    repo = Path(__file__).resolve().parents[1]
    wsl_repo = win_to_wsl(repo)
    inner = "cd " + shlex.quote(wsl_repo) + " && python3 -m mercury_merr.pyscf_run"
    if argv:
        inner += " " + " ".join(shlex.quote(a) for a in argv)
    cmd = ["wsl.exe", "-d", "Ubuntu-26.04", "--", "bash", "-lc", inner]
    print("reexec:", cmd, flush=True)
    return subprocess.call(cmd)


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(prog="pyscf-run")
    ap.add_argument("--structures", default="structures/species")
    ap.add_argument("--jobs", default="jobs/pyscf")
    ap.add_argument("--out", default="results/energies.csv")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--opt", action="store_true")
    args = ap.parse_args(argv)
    if not has_pyscf():
        extra = list(argv) if argv is not None else sys.argv[1:]
        return reexec_wsl(extra)
    os.environ.setdefault("OMP_NUM_THREADS", "8")
    run_all(args.structures, args.jobs, args.out, args.config, do_opt=args.opt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
