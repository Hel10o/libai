#!/usr/bin/env python3
"""Reproduce Q1 from the original input in a separate output directory.
Usage: python run_all.py --out reproduced
The complete default run includes all finite-volume, independent-reference,
conservation, sensitivity and workbook audits, and writes figures/paper.
"""
from __future__ import annotations
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
# Set before NumPy/SciPy imports (also inherited by all child processes).
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    b=Path(__file__).resolve().parent;ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=b/'input'/'附件1.xlsx')
    ap.add_argument('--template',type=Path,default=b/'input'/'result1_template.xlsx')
    ap.add_argument('--out',type=Path,default=b/'reproduced')
    args=ap.parse_args();args.out=args.out.resolve();args.out.mkdir(parents=True,exist_ok=True)
    if args.out.resolve()==(b/'input').resolve():raise ValueError('Output must not be the original-input folder')
    originals={str(p.resolve()):sha(p) for p in (args.input,args.template)}
    start=time.perf_counter();steps=[]
    for module,params in [
      ('q1_solver.py',['--input',str(args.input),'--out',str(args.out),'--quadrature']),
      ('q1_validate.py',['--input',str(args.input),'--out',str(args.out/'validation')]),
      ('q1_refine.py',['--input',str(args.input),'--out',str(args.out/'validation')]),
      ('q1_excel.py',['--template',str(args.template),'--raw',str(args.out/'q1_unrounded.npz'),'--out',str(args.out/'result1.xlsx')]),
      ('q1_report.py',['--input',str(args.input),'--out',str(args.out)])]:
        cmd=[sys.executable,str(b/'source'/module)]+params
        print('RUN',module,flush=True);tic=time.perf_counter()
        with (args.out/(module.replace('.py','')+'.log')).open('w',encoding='utf-8') as log:
            p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=os.environ.copy())
        steps.append({'module':module,'returncode':p.returncode,'elapsed_seconds':time.perf_counter()-tic})
        if p.returncode:raise RuntimeError(f'{module} failed; see {args.out/module.replace(".py",".log")}')
    for path,h in originals.items():assert sha(path)==h,'Original input changed'
    # Core data reproducibility compared with delivered numerical outputs if available.
    import numpy as np
    report={'all_steps_successful':True,'steps':steps,'elapsed_seconds':time.perf_counter()-start,
            'original_input_hashes_unchanged':True,'input_sha256':originals}
    delivered=b/'output'/'q1_unrounded.npz'
    if delivered.exists() and delivered.resolve()!=(args.out/'q1_unrounded.npz').resolve():
        a=np.load(delivered);c=np.load(args.out/'q1_unrounded.npz')
        report['compared_with_delivered']={k:{'max_abs_difference':float(np.max(np.abs(a[k]-c[k]))),
             'four_decimal_equal':bool(np.array_equal(np.round(a[k],4),np.round(c[k],4)))}
             for k in ('temperature_degC','moisture_dry_basis')}
    (args.out/'reproduction_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
if __name__=='__main__':main()
