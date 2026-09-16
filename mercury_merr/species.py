from __future__ import annotations
import math
import numpy as np
from .geometry import write_xyz

SPECIES_META={"hgcl2":(0,1),"hgcl3":(-1,1),"hgcl4":(-2,1),"hgoh2":(0,1),"chloride":(-1,1),"hydroxide":(-1,1),"methanethiolate":(-1,1),"hg_sme3":(-1,1)}

def _tetra_vectors():
    a=1.0/math.sqrt(3.0); return [np.array([a,a,a]),np.array([a,-a,-a]),np.array([-a,a,-a]),np.array([-a,-a,a])]

def build_species(name):
    if name=="hgcl2":
        d=2.35; return [("Hg",0,0,0),("Cl",0,0,d),("Cl",0,0,-d)]
    if name=="hgcl3":
        d=2.50; atoms=[("Hg",0,0,0)]
        for k in range(3):
            p=2*math.pi*k/3; atoms.append(("Cl",d*math.cos(p),d*math.sin(p),0))
        return atoms
    if name=="hgcl4":
        d=2.62; atoms=[("Hg",0,0,0)]
        for v in _tetra_vectors():
            r=d*v; atoms.append(("Cl",*r))
        return atoms
    if name=="hgoh2":
        d_hgo=2.10; d_oh=0.97; return [("Hg",0,0,0),("O",0,0,d_hgo),("H",0,0,d_hgo+d_oh),("O",0,0,-d_hgo),("H",0,0,-d_hgo-d_oh)]
    if name=="chloride": return [("Cl",0,0,0)]
    if name=="hydroxide": return [("O",0,0,0),("H",0,0,0.97)]
    if name=="methanethiolate":
        atoms=[("S",0,0,0),("C",1.82,0,0)]; c=np.array([1.82,0,0]); s=np.array([0,0,0]); z=(s-c)/np.linalg.norm(s-c); x=np.array([0,1,0.0]); y=np.cross(z,x); y/=np.linalg.norm(y); theta=math.radians(109.47); ch=1.09
        for k in range(3):
            phi=2*math.pi*k/3; v=math.cos(theta)*z+math.sin(theta)*(math.cos(phi)*x+math.sin(phi)*y); r=c+ch*v; atoms.append(("H",*r))
        return atoms
    if name=="hg_sme3":
        hg_s=2.45; s_c=1.82; atoms=[("Hg",0,0,0)]
        for k in range(3):
            p=2*math.pi*k/3; s=np.array([hg_s*math.cos(p),hg_s*math.sin(p),0]); outward=s/np.linalg.norm(s); c=s+s_c*outward; atoms.append(("S",*s)); atoms.append(("C",*c)); z=(s-c)/np.linalg.norm(s-c); trial=np.array([0,0,1.0]); x=np.cross(z,trial)
            if np.linalg.norm(x)<1e-8: x=np.array([1.0,0,0])
            x/=np.linalg.norm(x); y=np.cross(z,x); y/=np.linalg.norm(y); theta=math.radians(109.47); ch=1.09
            for j in range(3):
                phi=2*math.pi*j/3; v=math.cos(theta)*z+math.sin(theta)*(math.cos(phi)*x+math.sin(phi)*y); h=c+ch*v; atoms.append(("H",*h))
        return atoms
    raise KeyError(name)

def build_all(outdir):
    from pathlib import Path
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    for name,(charge,mult) in SPECIES_META.items(): write_xyz(outdir/f"{name}.xyz",build_species(name),f"{name}; charge={charge}; multiplicity={mult}; starting geometry")
