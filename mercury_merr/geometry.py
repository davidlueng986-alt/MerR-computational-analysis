from __future__ import annotations
import math
from pathlib import Path
import numpy as np

Atom = tuple[str, float, float, float]

def write_xyz(path: str | Path, atoms: list[Atom], comment: str = "") -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for el, x, y, z in atoms: f.write(f"{el:2s} {x: .8f} {y: .8f} {z: .8f}\n")

def read_xyz(path: str | Path) -> list[Atom]:
    lines = Path(path).read_text().splitlines(); n = int(lines[0]); atoms=[]
    for line in lines[2:2+n]:
        p=line.split(); atoms.append((p[0],float(p[1]),float(p[2]),float(p[3])))
    return atoms

def unit(v):
    v=np.asarray(v,dtype=float); n=np.linalg.norm(v)
    if n<1e-12: raise ValueError("zero vector")
    return v/n

def orthonormal_frame(axis):
    z=unit(axis); trial=np.array([1.0,0.0,0.0]) if abs(z[0])<0.8 else np.array([0.0,1.0,0.0]); x=unit(np.cross(trial,z)); y=unit(np.cross(z,x)); return x,y,z

def methyl_cap_hydrogens(carbon,sulfur,ch=1.09,angle_deg=109.47):
    c=np.array(carbon,float); s=np.array(sulfur,float); z=unit(s-c); x,y,_=orthonormal_frame(z); theta=math.radians(angle_deg); zs=math.cos(theta); rs=math.sin(theta); hs=[]
    for k in range(3):
        phi=2*math.pi*k/3.0; v=zs*z+rs*(math.cos(phi)*x+math.sin(phi)*y); hs.append(c+ch*v)
    return hs
