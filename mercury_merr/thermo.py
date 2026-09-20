from __future__ import annotations
from pathlib import Path
import math,yaml
import pandas as pd
import matplotlib.pyplot as plt
R_KCAL=0.00198720425864083

def load_energies(csv):
    df=pd.read_csv(csv); return {r.species:float(r.G_composite_kcal) for _,r in df.dropna(subset=["G_composite_kcal"]).iterrows()}

def reaction_dg(stoich,energies):
    missing=[s for s in stoich if s not in energies]
    if missing: return None,missing
    return sum(float(nu)*energies[s] for s,nu in stoich.items()),[]

def conditional_dg(dg0,T,chloride_products=0,chloride_M=1.0,hydroxide_products=0,pH=7.0,pKw=14.0):
    q=1.0
    if chloride_products: q*=max(chloride_M,1e-30)**chloride_products
    if hydroxide_products:
        oh=10**(-(pKw-pH)); q*=max(oh,1e-30)**hydroxide_products
    return dg0+R_KCAL*T*math.log(q)

def analyze(energies_csv,reactions_yaml,config_yaml,outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True); E=load_energies(energies_csv); rx=yaml.safe_load(Path(reactions_yaml).read_text())["reactions"]; cfg=yaml.safe_load(Path(config_yaml).read_text()); T=float(cfg["chemistry"]["temperature_K"])
    rows=[]
    for name,r in rx.items():
        dg0,missing=reaction_dg(r["stoich"],E); rows.append({"reaction":name,"description":r.get("description",""),"dG0_kcal_mol":dg0,"missing":",".join(missing)})
    pd.DataFrame(rows).to_csv(outdir/"reaction_dG0.csv",index=False)
    curves=[]
    for name,r in rx.items():
        if r.get("chloride_products",0)<=0: continue
        dg0,miss=reaction_dg(r["stoich"],E)
        if miss: continue
        for cl in cfg["chemistry"]["chloride_M_grid"]: curves.append({"reaction":name,"chloride_M":cl,"dG_cond_kcal_mol":conditional_dg(dg0,T,chloride_products=r.get("chloride_products",0),chloride_M=cl)})
    cdf=pd.DataFrame(curves); cdf.to_csv(outdir/"chloride_conditional_dG.csv",index=False)
    if not cdf.empty:
        fig,ax=plt.subplots(figsize=(7,5))
        for name,g in cdf.groupby("reaction"): ax.plot(g["chloride_M"],g["dG_cond_kcal_mol"],marker="o",label=name)
        ax.set_xscale("log"); ax.set_xlabel("[Cl-] (M; activity≈concentration screening assumption)"); ax.set_ylabel("Conditional ΔG (kcal/mol)")
        ax.set_title("engine=pyscf screening (not experimental Kd)")
        ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(outdir/"chloride_conditional_dG.png",dpi=200); fig.savefig(outdir/"chloride_conditional_dG.svg"); plt.close(fig)
    phrows=[]
    for name,r in rx.items():
        if r.get("hydroxide_products",0)<=0: continue
        dg0,miss=reaction_dg(r["stoich"],E)
        if miss: continue
        for ph in cfg["chemistry"]["pH_grid"]: phrows.append({"reaction":name,"pH":ph,"dG_cond_kcal_mol":conditional_dg(dg0,T,hydroxide_products=r.get("hydroxide_products",0),pH=ph)})
    pdf=pd.DataFrame(phrows); pdf.to_csv(outdir/"pH_conditional_dG.csv",index=False)
    if not pdf.empty:
        fig,ax=plt.subplots(figsize=(7,5))
        for name,g in pdf.groupby("reaction"): ax.plot(g["pH"],g["dG_cond_kcal_mol"],marker="o",label=name)
        ax.set_xlabel("pH (ideal Kw screening approximation)"); ax.set_ylabel("Conditional ΔG (kcal/mol)")
        ax.set_title("engine=pyscf screening (not experimental Kd)")
        ax.legend(); fig.tight_layout(); fig.savefig(outdir/"pH_conditional_dG.png",dpi=200); fig.savefig(outdir/"pH_conditional_dG.svg"); plt.close(fig)
