"""Reproduce the executed Jacobian and dimension-reduction checks in a new folder."""
import argparse,json
from pathlib import Path
import numpy as np
from geometry_solver import GeometryFV,Inputs

def checks(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    names=['jacobian_flux_checks.json','dimension_reduction_check.json']
    if any((out/n).exists() for n in names):raise FileExistsError('Use a new validation output directory')
    inp=Inputs();records={}
    for flux in ('integral','harmonic'):
        for case in ('q1','q23','q4'):
            op=GeometryFV(8,6,case,inp,flux=flux);rng=np.random.default_rng(3402)
            y=np.column_stack([28+10*rng.random(op.m),1.5+.5*rng.random(op.m)]).ravel()
            v=rng.standard_normal(2*op.m);eps=1e-5
            analytic=op.rhs(4000,y,inp.ambient,True)@v
            numeric=(op.rhs(4000,y+eps*v,inp.ambient)-op.rhs(4000,y-eps*v,inp.ambient))/(2*eps)
            r={'jacobian_direction_max_abs':float(np.max(abs(analytic-numeric))),
               'jacobian_direction_relative_norm':float(np.linalg.norm(analytic-numeric)/np.linalg.norm(analytic))}
            assert r['jacobian_direction_relative_norm']<1e-8
            records[case+'_'+flux]=r
    (out/names[0]).write_text(json.dumps(records,indent=2)+'\n')
    op=GeometryFV(20,16,'q4',inp,end_transfer=False);one=GeometryFV(20,0,'q4',inp)
    r=one.r;profile=np.column_stack([28+8*(r/.02)**2,2.55-.5*(r/.02)**2])
    y=np.broadcast_to(profile,(17,21,2)).copy().ravel()
    a=op.rhs(9000,y,inp.ambient).reshape(17,21,2)
    b=one.rhs(9000,profile.ravel(),inp.ambient).reshape(21,2)
    err=float(np.max(abs(a-b)));assert err<1e-12
    (out/names[1]).write_text(json.dumps({'no_end_exchange_radially_identical_rhs_max_abs':err,'time_s':9000,'case':'q4 moving material reference'},indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out-dir',required=True,type=Path)
    checks(p.parse_args().out_dir)
