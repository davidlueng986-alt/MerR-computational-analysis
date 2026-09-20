from __future__ import annotations

import json
from pathlib import Path
import re
import pandas as pd

from .species import SPECIES_META

HARTREE_TO_KCAL = 627.509474
PATTERNS = {
    "E_elec": re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"),
    "G_final": re.compile(r"Final Gibbs free energy\s+.*?(-?\d+\.\d+)\s+Eh", re.I),
    "G_corr": re.compile(r"G-E\(el\)\s+.*?(-?\d+\.\d+)\s+Eh", re.I),
}
ORCA_OK = "ORCA TERMINATED NORMALLY"
ORCA_METHOD = "ORCA PBE0-D4 ZORA (optfreq TZVP + SP TZVPP if present; not PySCF PBE0 DKH2 def2-TZVP)"


def last_match(pattern, text):
    m = pattern.findall(text)
    return float(m[-1]) if m else None


def _infer_engine(data: dict | None, method: str | None = None) -> str | None:
    if data:
        if data.get("engine"):
            return str(data["engine"])
        method = data.get("method", method)
    ml = str(method or "").lower()
    if "pyscf" in ml:
        return "pyscf"
    if "orca" in ml:
        return "orca"
    return None


def _empty_row(species: str) -> dict:
    return {
        "species": species,
        "E_opt_Eh": None,
        "G_opt_Eh": None,
        "Gcorr_Eh": None,
        "E_sp_Eh": None,
        "E_elec_Eh": None,
        "G_composite_Eh": None,
        "G_composite_kcal": None,
        "converged": False,
        "missing_sp": True,
        "g_source": None,
        "engine": None,
        "method": None,
        "error": None,
        "basis": None,
    }


def _from_energy_json(jobdir: Path, data: dict) -> dict:
    row = _empty_row(data.get("species", jobdir.name))
    e_elec = data.get("E_elec_Eh")
    g_comp = data.get("G_composite_Eh")
    g_kcal = data.get("G_composite_kcal")
    converged = bool(data.get("converged"))
    if not converged:
        g_comp = None
        g_kcal = None
    method = data.get("method")
    row.update(
        {
            "E_opt_Eh": data.get("E_opt_Eh"),
            "G_opt_Eh": data.get("G_opt_Eh"),
            "Gcorr_Eh": data.get("Gcorr_Eh"),
            "E_sp_Eh": data.get("E_sp_Eh", e_elec),
            "E_elec_Eh": e_elec,
            "G_composite_Eh": g_comp,
            "G_composite_kcal": g_kcal,
            "converged": converged,
            "missing_sp": bool(data.get("missing_sp", False)),
            "g_source": data.get("g_source", "energy.json"),
            "engine": _infer_engine(data, method),
            "method": method,
            "error": data.get("error"),
            "basis": data.get("basis"),
        }
    )
    return row


def parse_one(jobdir):
    jobdir = Path(jobdir)
    ej = jobdir / "energy.json"
    if ej.exists():
        try:
            data = json.loads(ej.read_text(encoding="utf-8"))
        except Exception as exc:
            row = _empty_row(jobdir.name)
            row["error"] = f"energy.json unreadable: {exc}"
            return row
        return _from_energy_json(jobdir, data)

    of = jobdir / "optfreq.out"
    if not of.exists():
        return None
    ot = of.read_text(errors="ignore")
    sp = jobdir / "sp.out"
    st = sp.read_text(errors="ignore") if sp.exists() else ""
    e_opt = last_match(PATTERNS["E_elec"], ot)
    g_final = last_match(PATTERNS["G_final"], ot)
    g_corr = last_match(PATTERNS["G_corr"], ot)
    e_sp = last_match(PATTERNS["E_elec"], st) if st else None
    if g_corr is None and g_final is not None and e_opt is not None:
        g_corr = g_final - e_opt

    opt_ok = ORCA_OK in ot
    sp_ok = bool(st) and ORCA_OK in st
    missing_sp = not sp_ok
    error = None
    if not opt_ok:
        error = "optfreq did not report ORCA TERMINATED NORMALLY"
    elif missing_sp:
        error = "missing or unconverged SP (TZVPP); G is optfreq TZVP if present"

    # Do not silently treat TZVP G as a TZVPP composite SP energy.
    if sp_ok and e_sp is not None and g_corr is not None:
        g_comp = e_sp + g_corr
        g_source = "sp_TZVPP+gcorr"
        e_elec = e_sp
    elif opt_ok and g_final is not None:
        g_comp = g_final
        g_source = "optfreq_TZVP_G"
        e_elec = e_opt
    else:
        g_comp = None
        g_source = None
        e_elec = e_sp if e_sp is not None else e_opt

    converged = bool(opt_ok and g_comp is not None)
    if not converged:
        g_comp = None

    return {
        "species": jobdir.name,
        "E_opt_Eh": e_opt,
        "G_opt_Eh": g_final,
        "Gcorr_Eh": g_corr,
        "E_sp_Eh": e_sp,
        "E_elec_Eh": e_elec,
        "G_composite_Eh": g_comp,
        "G_composite_kcal": None if g_comp is None else g_comp * HARTREE_TO_KCAL,
        "converged": converged,
        "missing_sp": missing_sp,
        "g_source": g_source,
        "engine": "orca",
        "method": ORCA_METHOD,
        "error": error,
        "basis": "ma-ZORA-def2-TZVP/TZVPP + SARC-ZORA",
    }


def parse_jobs(jobs_dir, out_csv):
    jobs_dir = Path(jobs_dir)
    rows = []
    if jobs_dir.exists():
        for d in sorted(jobs_dir.iterdir()):
            if not d.is_dir():
                continue
            if d.name in SPECIES_META:
                r = parse_one(d)
                if r:
                    rows.append(r)
                continue
            for sub in sorted(d.iterdir()):
                if sub.is_dir() and sub.name in SPECIES_META:
                    r = parse_one(sub)
                    if r:
                        rows.append(r)
    df = pd.DataFrame(rows)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df
