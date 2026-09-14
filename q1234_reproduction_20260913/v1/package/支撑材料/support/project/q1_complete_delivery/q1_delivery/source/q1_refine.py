#!/usr/bin/env python3
"""Further refinement and genuinely independent nonlinear collocation diagnostic."""
from pathlib import Path
import argparse, json, numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import BarycentricInterpolator
from scipy.optimize import brentq
from q1_solver import *
from q1_validate import differences, finite_cylinder_heat, TABLE_TIMES, ramp_cylinder

def chebyshev_moisture(env,N=128,rtol=3e-11,atol=3e-13):
    """Independent strong-form spectral differentiation in x=(r/R)^2.
    C_t=4/R² [D(C)*C_x+x*d_x(D(C)*C_x)]. Surface C is algebraic.
    Gauss-Lobatto nodes cluster near both axis and surface.
    Prescribed uniform initial state is retained at t=0 in outputs; the
    incompatible algebraic surface trace at 0+ is not substituted into it.
    """
    x=(1-np.cos(np.pi*np.arange(N+1)/N))/2
    w=(-1.)**np.arange(N+1);w[[0,-1]]*=.5
    dx=x[:,None]-x[None,:];np.fill_diagonal(dx,1.)
    mat=(w[None,:]/w[:,None])/dx;np.fill_diagonal(mat,0.)
    np.fill_diagonal(mat,-mat.sum(axis=1))
    assert np.max(np.abs(mat@x-1))<1e-9
    lastrow=mat[-1,:-1];dnn=mat[-1,-1]
    def full(t,y):
        cinf=env(t)[1];v=lastrow@y
        def g(c):return 2*P.D0*np.exp(-P.a/c)/P.R*(v+dnn*c)+P.hm*(c-cinf)
        cs=brentq(g,max(cinf+1e-8,0.5*float(np.min(y))),max(3.,2*cinf,-2*v/dnn),xtol=5e-15)
        return np.r_[y,cs]
    def rhs(t,y):
        c=full(t,y);cx=mat@c;flux=P.D0*np.exp(-P.a/c)*cx
        return (4/P.R**2*(flux+x*(mat@flux)))[:-1]
    y=np.full(N,P.C0);ts=np.arange(1801.);out=np.empty((1801,21));out[0]=P.C0
    xo=np.linspace(0,1,21)**2
    nf=0;steps=0
    for a,b in zip(np.arange(0,1800,60.),np.arange(60,1801,60.)):
        sol=solve_ivp(rhs,(a,b),y,method='BDF',rtol=rtol,atol=atol,max_step=6.,dense_output=True,first_step=min(.0001,.0001*(96/N)**4))
        if not sol.success: raise RuntimeError(sol.message)
        for t in np.arange(a+1,b+1):
            c=full(t,sol.sol(t));out[int(t)]=BarycentricInterpolator(x,c,wi=w)(xo)
        y=sol.y[:,-1];nf+=sol.nfev;steps+=len(sol.t)-1
    return out,{'N':N,'rtol':rtol,'atol':atol,'accepted_steps':steps,'nfev':nf,'method':'Independent Chebyshev collocation, algebraic surface Robin, BDF'}

def main():
    base=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=base/'input'/'附件1.xlsx')
    ap.add_argument('--out',type=Path,default=base/'output'/'validation')
    args=ap.parse_args();out=args.out
    data,_=read_xlsx(args.input);env=Environment(data)
    report=json.loads((out/'validation.json').read_text());R={}
    for n in [5120,10240]:
        R[n]={}
        for f in ['T','C']:
            if (out/f'{f}_tight_n{n}.npz').exists():
                z=dict(np.load(out/f'{f}_tight_n{n}.npz'));z['stats']=json.loads((out/f'{f}_tight_n{n}.npz.json').read_text())
            else:
                z=solve_fv(env,f,n=n,rtol=2e-12,atol=2e-14,max_step=3.,quadrature=True)
                save_field(out/f'{f}_tight_n{n}.npz',z)
            R[n][f]=z
        print('refined',n,flush=True)
    report['further_refinement']=[{'n_coarse':a,'n_fine':b,**{f:differences(R[a][f],R[b][f]) for f in ['T','C']}} for a,b in [(5120,10240)]]
    # Production uses 10240; compare 5120/10240 and independent spectral references.
    report['production_recommendation']={'n':10240,'rtol':2e-12,'atol':2e-14,'max_step_s':3.}
    tt={}
    for f in ['T','C']:
        z=solve_fv(env,f,n=10240,rtol=5e-13,atol=5e-15,max_step=1.)
        save_field(out/f'{f}_ultratight_n10240.npz',z)
        tt[f]=differences(R[10240][f],z)
    report['production_time_precision']=tt
    # Collocation cannot be mistaken for repeating the finite-volume routine.
    coll={}
    for N in [96,160,240]:
        coll[N],st=chebyshev_moisture(env,N)
        np.savez_compressed(out/f'independent_cheb_N{N}.npz',C=coll[N])
        print('cheb',N,differences(R[10240]['C'],coll[N]),flush=True)
    report['independent_nonlinear_collocation']={'N96_vs_N160':differences(coll[96],coll[160]),
      'N160_vs_N240':differences(coll[160],coll[240]),'fv10240_vs_N240':differences(R[10240]['C'],coll[240]),
      'settings':st}
    a,ma=finite_cylinder_heat(data,TABLE_TIMES,R[10240]['T']['radius_m'],[0,P.L/2-.02,P.L/2],200,200)
    b,mb=finite_cylinder_heat(data,TABLE_TIMES,R[10240]['T']['radius_m'],[0,P.L/2-.02,P.L/2],400,400)
    report['2d_heat_refinement']={'200_vs_400_midplane_max_K':float(np.max(np.abs(a[:,:,0]-b[:,:,0]))),
       '200_vs_400_all_max_K':float(np.max(np.abs(a-b))),
       '400_midplane_vs_production_K':float(np.max(np.abs(b[:,:,0]-R[10240]['T']['values'][TABLE_TIMES]))),
       'average_2d_1800':float(mb[-1]),'end_center_1800':float(b[-1,0,-1]),'end_corner_1800':float(b[-1,-1,-1])}
    np.savez_compressed(out/'finite_cylinder_2d_heat_refined.npz',time_s=TABLE_TIMES,radius_m=R[10240]['T']['radius_m'],z_m=[0,P.L/2-.02,P.L/2],T=b,average=mb)
    thermal=ramp_cylinder(data,R[10240]['T']['time_s'],R[10240]['T']['radius_m'],800)
    report['production_actual_thermal_analytic']=differences(R[10240]['T'],thermal)
    (out/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['further_refinement','production_time_precision','independent_nonlinear_collocation','2d_heat_refinement','production_actual_thermal_analytic']},indent=2),flush=True)
if __name__=='__main__':main()
