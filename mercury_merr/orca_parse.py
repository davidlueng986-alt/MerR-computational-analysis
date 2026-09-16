from __future__ import annotations
from pathlib import Path
import re
import pandas as pd

HARTREE_TO_KCAL=627.509474
PATTERNS={"E_elec":re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"),"G_final":re.compile(r"Final Gibbs free energy\s+.*?(-?\d+\.\d+)\s+Eh",re.I),"G_corr":re.compile(r"G-E\(el\)\s+.*?(-?\d+\.\d+)\s+Eh",re.I)}

def last_match(pattern,text):
    m=pattern.findall(text); return float(m[-1]) if m else None

def parse_one(jobdir):
    jobdir=Path(jobdir); of=jobdir/"optfreq.out"; sp=jobdir/"sp.out"
    if not of.exists(): return None
    ot=of.read_text(errors="ignore"); st=sp.read_text(errors="ignore") if sp.exists() else ""
    e_opt=last_match(PATTERNS["E_elec"],ot); g_final=last_match(PATTERNS["G_final"],ot); g_corr=last_match(PATTERNS["G_corr"],ot); e_sp=last_match(PATTERNS["E_elec"],st) if st else None
    if g_corr is None and g_final is not None and e_opt is not None: g_corr=g_final-e_opt
    g_comp=(e_sp+g_corr) if (e_sp is not None and g_corr is not None) else g_final
    return {"species":jobdir.name,"E_opt_Eh":e_opt,"G_opt_Eh":g_final,"Gcorr_Eh":g_corr,"E_sp_Eh":e_sp,"G_composite_Eh":g_comp,"G_composite_kcal":None if g_comp is None else g_comp*HARTREE_TO_KCAL}

def parse_jobs(jobs_dir,out_csv):
    rows=[]
    for d in sorted(Path(jobs_dir).iterdir()):
        if d.is_dir():
            r=parse_one(d)
            if r: rows.append(r)
    df=pd.DataFrame(rows); df.to_csv(out_csv,index=False); return df
