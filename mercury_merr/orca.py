from __future__ import annotations
from pathlib import Path
import yaml
from .species import SPECIES_META


def load_config(path): return yaml.safe_load(Path(path).read_text())


def orca_input(xyz_path,charge,mult,keywords,hg_basis,nprocs,maxcore):
    return f'''! {keywords} PrintBasis\n\n%pal\n  nprocs {nprocs}\nend\n%maxcore {maxcore}\n%basis\n  NewGTO Hg "{hg_basis}" end\n  DelECP Hg\nend\n\n* xyzfile {charge} {mult} {xyz_path}\n'''


def generate_jobs(config_path,structures_dir,jobs_dir):
    cfg=load_config(config_path); oc=cfg["orca"]; structures_dir=Path(structures_dir).resolve(); jobs_dir=Path(jobs_dir); jobs_dir.mkdir(parents=True,exist_ok=True)
    for name,(charge,mult) in SPECIES_META.items():
        xyz=structures_dir/f"{name}.xyz"; d=jobs_dir/name; d.mkdir(parents=True,exist_ok=True)
        (d/"optfreq.inp").write_text(orca_input(str(xyz),charge,mult,oc["optfreq_keywords"],oc["hg_basis_opt"],oc["nprocs"],oc["maxcore_mb"]))
        (d/"sp.inp").write_text(orca_input("optfreq.xyz",charge,mult,oc["sp_keywords"],oc["hg_basis_sp"],oc["nprocs"],oc["maxcore_mb"]))
        (d/"run.sh").write_text(f'''#!/usr/bin/env bash\nset -euo pipefail\nORCA=${{ORCA:-{oc["executable"]}}}\n"$ORCA" optfreq.inp > optfreq.out\ngrep -q "ORCA TERMINATED NORMALLY" optfreq.out\n"$ORCA" sp.inp > sp.out\ngrep -q "ORCA TERMINATED NORMALLY" sp.out\n'''); (d/"run.sh").chmod(0o755)
