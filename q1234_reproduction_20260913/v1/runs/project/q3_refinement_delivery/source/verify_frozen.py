"""Read-only verifier for a frozen q3_refinement_delivery directory."""
from __future__ import annotations
import argparse, hashlib, json, os, sys
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
from validate_xlsx import validate as validate_xlsx


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def snapshot(root: Path):
    return {str(p.relative_to(root)):sha256(p) for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}


def read_manifest(path:Path):
    out={}
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.strip():continue
        digest,rel=line.split('  ',1);out[rel]=digest
    return out


def verify(root:Path,audit_out:Path):
    root=root.resolve();audit_out=audit_out.resolve()
    try:audit_out.relative_to(root)
    except ValueError:pass
    else:raise ValueError('audit-out must be outside the frozen root')
    if audit_out.exists():raise FileExistsError(audit_out)
    before=snapshot(root)
    manifest=read_manifest(root/'MANIFEST.sha256')
    actual={k:v for k,v in before.items() if k!='MANIFEST.sha256'}
    if set(manifest)!=set(actual):
        raise ValueError(f'manifest coverage mismatch: missing={sorted(set(actual)-set(manifest))[:5]}, stale={sorted(set(manifest)-set(actual))[:5]}')
    bad=[k for k in manifest if manifest[k]!=actual[k]]
    if bad:raise ValueError(f'manifest mismatch: {bad[:8]}')
    xlsx=validate_xlsx(root/'output/result3.xlsx',root/'output/solution.npz')
    event=json.loads((root/'output/end_event.json').read_text(encoding='utf-8'))
    z=np.load(root/'output/solution.npz')
    maxC=float(z['profile_TC'][-1,:,1].max())
    checks={
      'event_execution_matches_solution':abs(float(z['time_s'][-1])-event['execution_s'])<1e-8,
      'strict_endpoint':maxC<0.15,
      'event_execution_max_matches_profile':abs(maxC-event['execution_max'])<1e-12,
      'table5_exists':(root/'output/table5.csv').is_file(),
      'xlsx':xlsx,
    }
    if not all(v if isinstance(v,bool) else True for v in checks.values()):raise AssertionError(checks)
    after=snapshot(root)
    if before!=after:raise AssertionError('Frozen root changed during verification')
    audit_out.mkdir(parents=True)
    obj={'passed':True,'root':str(root),'root_unchanged':True,'manifest_files':len(manifest),'endpoint_max_C':maxC,'checks':checks}
    (audit_out/'verify.json').write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return obj

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--audit-out',type=Path,required=True);a=p.parse_args();print(json.dumps(verify(a.root,a.audit_out),ensure_ascii=False,indent=2))
