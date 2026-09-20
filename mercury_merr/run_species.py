"""CLI: python -m mercury_merr.run_species --name <species> --engine pyscf --workdir ."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from .species import SPECIES_META, build_all


def win_to_wsl(path: Path) -> str:
    p = str(path.resolve())
    if len(p) >= 2 and p[1] == ":":
        return "/mnt/" + p[0].lower() + p[2:].replace("\\", "/")
    return p.replace("\\", "/")


def reexec_wsl(argv: list[str], workdir: Path) -> int:
    wsl_wd = win_to_wsl(workdir)
    inner = "cd " + shlex.quote(wsl_wd) + " && python3 -m mercury_merr.run_species"
    if argv:
        inner += " " + " ".join(shlex.quote(a) for a in argv)
    distros = ["Ubuntu-26.04", "Ubuntu", "Debian"]
    last_err = None
    for distro in distros:
        cmd = ["wsl.exe", "-d", distro, "--", "bash", "-lc", inner]
        print("reexec:", cmd, flush=True)
        try:
            return subprocess.call(cmd)
        except Exception as exc:
            last_err = exc
    raise RuntimeError(f"WSL reexec failed: {last_err}")


def _write_failed(jobdir: Path, name: str, engine: str, error: str) -> dict:
    from .pyscf_engine import _blank, write_energy_json

    jobdir.mkdir(parents=True, exist_ok=True)
    rec = _blank(name, *SPECIES_META[name], "n/a")
    rec["method"] = engine
    rec["error"] = error
    rec["converged"] = False
    write_energy_json(jobdir, rec)
    return rec


def run_pyscf_species(name: str, workdir: Path, reuse_opt: bool = False) -> dict:
    from .geometry import write_xyz
    from .pyscf_engine import run_species
    from .species import build_species

    structures = workdir / "structures" / "species"
    structures.mkdir(parents=True, exist_ok=True)
    xyz = structures / f"{name}.xyz"
    charge, mult = SPECIES_META[name]
    write_xyz(xyz, build_species(name), f"{name}; charge={charge}; multiplicity={mult}; starting geometry")
    return run_species(name, xyz, workdir / "jobs" / name, reuse_opt=reuse_opt)


def run_orca_species(name: str, workdir: Path) -> dict:
    from .orca import generate_jobs

    structures = workdir / "structures" / "species"
    if not (structures / f"{name}.xyz").exists():
        build_all(structures)
    generate_jobs(workdir / "config" / "config.yaml", structures, workdir / "jobs")
    jobdir = workdir / "jobs" / name
    exe = shutil.which("orca") or shutil.which("orca.exe")
    if not exe:
        return _write_failed(
            jobdir,
            name,
            "ORCA PBE0-D4 ZORA",
            "ORCA executable not found; use --engine pyscf",
        )
    run_sh = jobdir / "run.sh"
    if not run_sh.exists():
        return _write_failed(jobdir, name, "ORCA PBE0-D4 ZORA", "missing run.sh")
    rc = subprocess.call(["bash", str(run_sh)], cwd=str(jobdir))
    if rc != 0:
        return _write_failed(jobdir, name, "ORCA PBE0-D4 ZORA", f"run.sh exited {rc}")
    from .orca_parse import parse_one

    row = parse_one(jobdir)
    (jobdir / "energy.json").write_text(json.dumps(row, indent=2, default=str) + "\n")
    return row


def _finalize(workdir: Path) -> None:
    from .orca_parse import parse_jobs
    from .thermo import analyze

    energies = workdir / "results" / "energies.csv"
    parse_jobs(workdir / "jobs", energies)
    rx = workdir / "config" / "reactions.yaml"
    cfg = workdir / "config" / "config.yaml"
    if rx.exists() and cfg.exists() and energies.exists():
        analyze(energies, rx, cfg, workdir / "results" / "thermo")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="run-species")
    ap.add_argument("--name", required=True, help="species name or 'all'")
    ap.add_argument("--engine", default="pyscf", choices=["pyscf", "orca"])
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--reuse-opt", action="store_true", help="reuse jobs/<name>/opt.xyz if present")
    args = ap.parse_args(argv)
    workdir = Path(args.workdir).resolve()
    names = list(SPECIES_META) if args.name in ("all", "*") else [args.name]
    for n in names:
        if n not in SPECIES_META:
            raise SystemExit(f"unknown species {n}; known: {', '.join(SPECIES_META)}")

    if args.engine == "pyscf":
        from .pyscf_engine import has_pyscf

        if not has_pyscf():
            extra = list(argv) if argv is not None else sys.argv[1:]
            return reexec_wsl(extra, workdir)
        os.environ.setdefault("OMP_NUM_THREADS", os.environ.get("OMP_NUM_THREADS", "8"))

    rows = []
    for n in names:
        print(f"=== {args.engine} {n} ===", flush=True)
        if args.engine == "pyscf":
            rec = run_pyscf_species(n, workdir, reuse_opt=args.reuse_opt)
        else:
            rec = run_orca_species(n, workdir)
        rows.append(rec)
        print(
            json.dumps(
                {k: rec.get(k) for k in ("species", "converged", "E_elec_Eh", "G_composite_kcal", "error", "seconds", "method")},
                default=str,
            ),
            flush=True,
        )
        try:
            _finalize(workdir)
        except Exception as exc:
            print(f"finalize warning after {n}: {exc}", flush=True)
    n_ok = sum(1 for r in rows if r.get("converged"))
    meta = {
        "engine": args.engine,
        "method": rows[0].get("method") if rows else None,
        "n_converged": n_ok,
        "n_species": len(rows),
        "energies": str(workdir / "results" / "energies.csv"),
    }
    (workdir / "results").mkdir(parents=True, exist_ok=True)
    (workdir / "results" / "qm_engine.json").write_text(json.dumps(meta, indent=2) + "\n")
    return 0 if n_ok == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
