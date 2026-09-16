from pathlib import Path
import requests

def fetch_pdb(pdb_id, out_path, template="https://files.rcsb.org/download/{pdb_id}.pdb"):
    url=template.format(pdb_id=pdb_id.upper())
    r=requests.get(url,timeout=60); r.raise_for_status()
    Path(out_path).write_bytes(r.content)
    return url
