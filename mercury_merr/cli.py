from __future__ import annotations
import argparse
from pathlib import Path
import yaml
from .fetch import fetch_pdb
from .species import build_all
from .pdbtools import donor_report, build_methylthiolate_site_from_pdb
from .orca import generate_jobs
from .orca_parse import parse_jobs
from .pyscf_run import main as pyscf_main
from .run_species import main as run_species_main
from .thermo import analyze


def main():
    ap=argparse.ArgumentParser(prog="merr-hg")
    sub=ap.add_subparsers(dest="cmd",required=True)
    p=sub.add_parser("fetch"); p.add_argument("--config",default="config/config.yaml"); p.add_argument("--out",default="structures/5CRL.pdb")
    p=sub.add_parser("build-species"); p.add_argument("--outdir",default="structures/species")
    p=sub.add_parser("inspect-site"); p.add_argument("pdb"); p.add_argument("--cutoff",type=float,default=3.2)
    p=sub.add_parser("build-site-mimic"); p.add_argument("pdb"); p.add_argument("--out",default="structures/5CRL_site_mimic.xyz"); p.add_argument("--site",type=int,default=0); p.add_argument("--cutoff",type=float,default=3.2)
    p=sub.add_parser("gen-orca"); p.add_argument("--config",default="config/config.yaml"); p.add_argument("--structures",default="structures/species"); p.add_argument("--jobs",default="jobs")
    p=sub.add_parser("parse-orca"); p.add_argument("--jobs",default="jobs"); p.add_argument("--out",default="results/energies.csv")
    p=sub.add_parser("run-pyscf"); p.add_argument("--structures",default="structures/species"); p.add_argument("--jobs",default="jobs/pyscf"); p.add_argument("--out",default="results/energies.csv"); p.add_argument("--config",default="config/config.yaml"); p.add_argument("--opt",action="store_true")
    p=sub.add_parser("run-species"); p.add_argument("--name",required=True); p.add_argument("--engine",default="pyscf",choices=["pyscf","orca"]); p.add_argument("--workdir",default="."); p.add_argument("--reuse-opt",action="store_true")
    p=sub.add_parser("thermo"); p.add_argument("--energies",default="results/energies.csv"); p.add_argument("--reactions",default="config/reactions.yaml"); p.add_argument("--config",default="config/config.yaml"); p.add_argument("--outdir",default="results/thermo")
    p=sub.add_parser("screen"); p.add_argument("--outdir",default="structures/species"); p.add_argument("--jobs",default="jobs/pyscf"); p.add_argument("--energies",default="results/energies.csv"); p.add_argument("--config",default="config/config.yaml"); p.add_argument("--reactions",default="config/reactions.yaml"); p.add_argument("--thermo-outdir",default="results/thermo"); p.add_argument("--opt",action="store_true")
    a=ap.parse_args()
    if a.cmd=="fetch":
        cfg=yaml.safe_load(Path(a.config).read_text()); pr=cfg["project"]; Path(a.out).parent.mkdir(parents=True,exist_ok=True); print(fetch_pdb(pr["pdb_id"],a.out,pr["pdb_url_template"]))
    elif a.cmd=="build-species": build_all(a.outdir)
    elif a.cmd=="inspect-site":
        for i,(hg,donors) in enumerate(donor_report(a.pdb,a.cutoff)):
            print(f"site {i}: Hg serial {hg['serial']} chain={hg['chain']} resseq={hg['resseq']}")
            for d,dd in donors: print(f"  {d['chain']}:{d['resname']}{d['resseq']} {d['name']} {dd:.3f} A")
    elif a.cmd=="build-site-mimic": print(build_methylthiolate_site_from_pdb(a.pdb,a.out,a.site,a.cutoff))
    elif a.cmd=="gen-orca": generate_jobs(a.config,a.structures,a.jobs)
    elif a.cmd=="parse-orca": Path(a.out).parent.mkdir(parents=True,exist_ok=True); print(parse_jobs(a.jobs,a.out))
    elif a.cmd=="run-pyscf":
        argv=["--structures",a.structures,"--jobs",a.jobs,"--out",a.out,"--config",a.config]
        if a.opt: argv.append("--opt")
        raise SystemExit(pyscf_main(argv))
    elif a.cmd=="run-species":
        argv=["--name",a.name,"--engine",a.engine,"--workdir",a.workdir]
        if a.reuse_opt: argv.append("--reuse-opt")
        raise SystemExit(run_species_main(argv))
    elif a.cmd=="thermo": analyze(a.energies,a.reactions,a.config,a.outdir)
    elif a.cmd=="screen":
        build_all(a.outdir)
        argv=["--structures",a.outdir,"--jobs",a.jobs,"--out",a.energies,"--config",a.config]
        if a.opt: argv.append("--opt")
        rc=pyscf_main(argv)
        if rc: raise SystemExit(rc)
        analyze(a.energies,a.reactions,a.config,a.thermo_outdir)

if __name__=="__main__": main()
