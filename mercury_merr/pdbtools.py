from __future__ import annotations
from pathlib import Path
import math
import numpy as np
from .geometry import write_xyz, methyl_cap_hydrogens

def parse_pdb(path):
    atoms=[]
    for line in Path(path).read_text().splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")): continue
        try:
            atom={"record":line[:6].strip(),"serial":int(line[6:11]),"name":line[12:16].strip(),"resname":line[17:20].strip(),"chain":line[21].strip(),"resseq":int(line[22:26]),"x":float(line[30:38]),"y":float(line[38:46]),"z":float(line[46:54]),"element":line[76:78].strip().upper() or line[12:16].strip()[0].upper(),"line":line}
        except Exception: continue
        atoms.append(atom)
    return atoms

def dist(a,b): return math.sqrt((a["x"]-b["x"])**2+(a["y"]-b["y"])**2+(a["z"]-b["z"])**2)
def find_hg_sites(atoms): return [a for a in atoms if a["element"]=="HG" or a["name"].upper()=="HG" or a["resname"].upper()=="HG"]

def donor_report(pdb_path,cutoff=3.2):
    atoms=parse_pdb(pdb_path); report=[]
    for hg in find_hg_sites(atoms):
        donors=[]
        for a in atoms:
            if a["resname"]=="CYS" and a["name"]=="SG":
                d=dist(hg,a)
                if d<=cutoff: donors.append((a,d))
        report.append((hg,sorted(donors,key=lambda x:x[1])))
    return report

def build_methylthiolate_site_from_pdb(pdb_path,out_xyz,site_index=0,cutoff=3.2):
    atoms=parse_pdb(pdb_path); reports=donor_report(pdb_path,cutoff)
    if not reports: raise RuntimeError("No Hg atoms found")
    hg,donors=reports[site_index]
    if len(donors)<3: raise RuntimeError(f"Expected >=3 Cys SG donors near Hg; found {len(donors)}. Increase cutoff or inspect structure.")
    donors=donors[:3]; out=[("Hg",hg["x"],hg["y"],hg["z"])]
    for sg,_ in donors:
        cb=next((a for a in atoms if a["chain"]==sg["chain"] and a["resseq"]==sg["resseq"] and a["name"]=="CB"),None)
        if cb is None: raise RuntimeError(f"CB not found for {sg['chain']}:{sg['resseq']}")
        s=np.array([sg["x"],sg["y"],sg["z"]]); c=np.array([cb["x"],cb["y"],cb["z"]]); out.append(("S",*s)); out.append(("C",*c))
        for h in methyl_cap_hydrogens(c,s): out.append(("H",*h))
    write_xyz(out_xyz,out,"5CRL-derived Hg(SMe)3 site mimic; optimize before interpretation")
    return [(d[0]["chain"],d[0]["resseq"],d[1]) for d in donors]
