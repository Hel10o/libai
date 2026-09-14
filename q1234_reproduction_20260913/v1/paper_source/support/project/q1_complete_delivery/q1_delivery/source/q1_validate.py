#!/usr/bin/env python3
"""Actual numerical tests, independent Bessel series, and model sensitivities."""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.special import j0,j1,jn_zeros
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from q1_solver import P,read_xlsx,Environment,ConstantEnvironment,RadialFV,solve_fv,save_field
TABLE_TIMES=np.array([100,300,600,900,1200,1500,1800])

def differences(a,b):
    x=a['values'] if isinstance(a,dict) else a
    y=b['values'] if isinstance(b,dict) else b
    d=np.abs(x-y)
    return {'max_full':float(d[1:].max()),'rms_full':float(np.sqrt(np.mean(d[1:]**2))),
            'max_table':float(d[TABLE_TIMES][:,::5].max()),
            'four_decimal_changed_full':int(np.sum(np.round(x[1:],4)!=np.round(y[1:],4))),
            'four_decimal_changed_table':int(np.sum(np.round(x[TABLE_TIMES][:,::5],4)!=np.round(y[TABLE_TIMES][:,::5],4))),
            'max_location_t_rindex':list(map(int,np.unravel_index(np.argmax(d),d.shape)))}

def roots_cylinder(Bi,modes):
    zj0=jn_zeros(0,modes);zj1=jn_zeros(1,modes)
    f=lambda z:z*j1(z)-Bi*j0(z)
    return np.array([brentq(f,1e-12 if i==0 else zj1[i-1]+1e-10,zj0[i]-1e-10,xtol=1e-13) for i in range(modes)])

def bessel_modes(diff,beta,modes=300):
    mu=roots_cylinder(beta*P.R/diff,modes)
    a=2*j1(mu)/(mu*(j0(mu)**2+j1(mu)**2))
    return mu,a,diff*(mu/P.R)**2

def constant_cylinder(times,r,initial,bath,diff,beta,modes=300):
    """Independently derived separation-of-variables Robin cylinder series."""
    mu,a,decay=bessel_modes(diff,beta,modes)
    phi=a[:,None]*j0(mu[:,None]*np.asarray(r)[None,:]/P.R)
    ans=bath+(initial-bath)*(np.exp(-np.asarray(times)[:,None]*decay)@phi)
    ans[np.asarray(times)==0]=initial
    return ans

def ramp_cylinder(data,times,r,modes=300):
    """Exact convolution on piecewise-linear environmental segments; no FV code."""
    mu,a,decay=bessel_modes(P.alpha,P.hT/(P.rho*P.cp),modes)
    phi=a[:,None]*j0(mu[:,None]*np.asarray(r)[None,:]/P.R)
    integ=np.zeros(modes);out=np.empty((len(times),len(r)));out[0]=P.T0
    env=Environment(data)
    for k,(ta,tb) in enumerate(zip(times[:-1],times[1:]),1):
        cuts=np.r_[ta,data[(data[:,0]>ta)&(data[:,0]<tb),0],tb]
        for l,h in zip(cuts[:-1],cuts[1:]):
            dt=h-l;slope=(env(h)[0]-env(l)[0])/dt
            ex=np.exp(-decay*dt)
            integ=integ*ex+slope*(-np.expm1(-decay*dt))/decay
        out[k]=env(tb)[0]-integ@phi
    return out

def finite_cylinder_heat(data,times,r,z,mr=160,mz=160):
    """2D axisymmetric r-z heat solution by product eigenfunction series.
    Half-domain z in [0,L/2]; same Robin hT on side and both ends.
    This is an actually evaluated finite-cylinder 2D solution, not 2D mass.
    """
    H=P.L/2;mu,A,dr=bessel_modes(P.alpha,P.hT/(P.rho*P.cp),mr)
    biz=P.hT*H/P.k
    lam=np.array([brentq(lambda x:x*np.tan(x)-biz,i*np.pi+1e-10,i*np.pi+np.pi/2-1e-10,xtol=1e-13) for i in range(mz)])
    B=4*np.sin(lam)/(2*lam+np.sin(2*lam))
    decay=dr[:,None]+P.alpha*(lam[None,:]/H)**2
    pr=A[:,None]*j0(mu[:,None]*np.asarray(r)[None,:]/P.R)
    pz=B[:,None]*np.cos(lam[:,None]*np.asarray(z)[None,:]/H)
    rav=A*2*j1(mu)/mu;zav=B*np.sin(lam)/lam
    I=np.zeros((mr,mz));out=[];means=[];env=Environment(data);last=0.
    for t in times:
        cuts=np.unique(np.r_[last,data[(data[:,0]>last)&(data[:,0]<t),0],t])
        for a,b in zip(cuts[:-1],cuts[1:]):
            dt=b-a;s=(env(b)[0]-env(a)[0])/dt
            I=I*np.exp(-decay*dt)+s*(-np.expm1(-decay*dt))/decay
        out.append(env(t)[0]-pr.T@I@pz)
        means.append(env(t)[0]-rav@I@zav);last=t
    return np.asarray(out),np.asarray(means)

class LatentEnvironment:
    """CONDITIONAL boundary diagnostic: all outward water loss evaporates at side.
    rho_d=initial_wet_density/(1+C0) is an EXTRA assumption, not appendix fact.
    lambda=2.45e6 J/kg is FAO56 Annex3 constant near 20degC.
    Energy: -k*T_r=hT*(Ts-Tinf)+lambda*rho_d*hm*(Cs-Cinf).
    """
    def __init__(self,env,time_s,Csurface,latent=2.45e6):
        self.env=env;self.cs=PchipInterpolator(time_s,Csurface);self.latent=latent
        self.rho_d=P.rho/(1+P.C0);self.kind='conditional_surface_latent'
        # Diagnostic only: boundary trace is interpolated at 1 s, not a startup-accurate latent reference.
        self.knots=np.unique(np.r_[env.knots, time_s[time_s<=10]])
    def __call__(self,t):
        a=self.env(t).copy();q=self.latent*self.rho_d*P.hm*(self.cs(t)-a[...,1])
        a[...,0]-=q/P.hT
        return a

def get_or_run(out,env,field,n,**kwargs):
    f=out/f'{field}_n{n}.npz'
    if f.exists() and not kwargs and 'first_step_policy' in json.loads(Path(str(f)+'.json').read_text()):
        z=dict(np.load(f));z['stats']=json.loads(Path(str(f)+'.json').read_text());return z
    z=solve_fv(env,field,n=n,**kwargs);save_field(f,z);return z

def main():
    ap=argparse.ArgumentParser();base=Path(__file__).resolve().parents[1]
    ap.add_argument('--input',type=Path,default=base/'input'/'附件1.xlsx')
    ap.add_argument('--out',type=Path,default=base/'output'/'validation')
    ap.add_argument('--grids',type=int,nargs='+',default=[20,40,80,160,320,640,1280,2560,5120])
    args=ap.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    data,audit=read_xlsx(args.input);env=Environment(data);report={'input_audit':audit,'tests':{},'grid':[]}
    cache={}
    for n in args.grids:
        cache[n]={}
        for f in ['T','C']:
            tic=time.perf_counter();cache[n][f]=get_or_run(out,env,f,n)
            print('grid',n,f,'elapsed',round(time.perf_counter()-tic,2),flush=True)
        report['grid'].append({'n':n,'dr_cm':100*P.R/n,**{f:cache[n][f]['stats'] for f in ['T','C']}})
    nf=args.grids[-1];ref=cache[nf]
    for row in report['grid']:
        row['difference_to_finest']={f:differences(cache[row['n']][f],ref[f]) for f in ['T','C']}
    report['grid_pairwise']=[{'n_coarse':a,'n_fine':b,**{f:differences(cache[a][f],cache[b][f]) for f in ['T','C']}}
                             for a,b in zip(args.grids[:-1],args.grids[1:])]
    (out/'progress.json').write_text(json.dumps(report,indent=2))
    # Time error isolated on a fixed finest spatial grid.
    tight={}
    for f in ['T','C']:
        tight[f]=solve_fv(env,f,n=nf,rtol=2e-12,atol=2e-14,max_step=3.,quadrature=True)
        save_field(out/f'{f}_tight_n{nf}.npz',tight[f])
        report['tests']['time_'+f]=differences(ref[f],tight[f])
        print('tight',f,report['tests']['time_'+f],flush=True)
    # Maximum-principle and actual boundary flux signs.
    for f in ['T','C']:
        v=tight[f]['values'];amb=env(tight[f]['time_s'])[:,0 if f=='T' else 1]
        report['tests']['physical_'+f]={'min':float(v.min()),'max':float(v.max()),
          'initial_error':float(np.max(np.abs(v[0]-(P.T0 if f=='T' else P.C0)))),
          'outward_boundary_difference_min':float(np.min(v[:,-1]-amb)),
          'outward_boundary_difference_max':float(np.max(v[:,-1]-amb)),
          'radial_monotonicity_violation':float(max(0,-np.min(np.diff(v,axis=1)))) if f=='T' else float(max(0,np.max(np.diff(v,axis=1))))}
    # No external driving; and arbitrary nonuniform profile, zero flux.
    for f in ['T','C']:
        ini=P.T0 if f=='T' else P.C0
        u=solve_fv(ConstantEnvironment(),f,n=80,end=600.,times=np.arange(601.))
        q=solve_fv(ConstantEnvironment(),f,n=80,end=600.,times=np.arange(601.),transfer_scale=0,
                   initial=lambda r:ini+0.2*np.cos(np.pi*r/P.R))
        report['tests']['uniform_'+f]={'max_error':float(np.max(np.abs(u['values']-ini))), 'pass':bool(np.max(np.abs(u['values']-ini))<1e-10)}
        err=float(np.max(np.abs(q['average']-q['average'][0])))
        report['tests']['no_flux_'+f]={'integral_average_drift':err,'pass':bool(err<1e-10)}
    # Constant-coefficient analytic check (new problem, same boundary implementation).
    diff0=P.D0*np.exp(-P.a/P.C0);ncheck=1280
    for f,diff,beta,ini,bath in [('T',P.alpha,P.hT/(P.rho*P.cp),P.T0,45.),('C',diff0,P.hm,P.C0,.025)]:
        ce=ConstantEnvironment(T=45,C=.025)
        u=solve_fv(ce,f,n=ncheck,constant_D=diff0 if f=='C' else None)
        a=constant_cylinder(u['time_s'],u['radius_m'],ini,bath,diff,beta,modes=400)
        b=constant_cylinder(u['time_s'],u['radius_m'],ini,bath,diff,beta,modes=800)
        report['tests']['constant_analytic_'+f]={'fv_vs_series':differences(u,b),'series_truncation':differences(a,b)}
        np.savez_compressed(out/f'constant_analytic_{f}.npz',time_s=u['time_s'],radius_m=u['radius_m'],fv=u['values'],series=b)
    # Independent direct solution of the actual thermal task.
    a=ramp_cylinder(data,tight['T']['time_s'],tight['T']['radius_m'],modes=300)
    b=ramp_cylinder(data,tight['T']['time_s'],tight['T']['radius_m'],modes=600)
    report['tests']['actual_thermal_analytic']={'fv_vs_series':differences(tight['T'],b),'series_truncation':differences(a,b)}
    np.savez_compressed(out/'actual_thermal_analytic.npz',temperature=b)
    print('analytical done',flush=True)
    # Finite-cylinder two-dimensional thermal comparison, no unperformed 2D claim.
    r=tight['T']['radius_m'];z=np.array([0.,P.L/2-.02,P.L/2])
    a,ma=finite_cylinder_heat(data,TABLE_TIMES,r,z,mr=100,mz=100)
    b,mb=finite_cylinder_heat(data,TABLE_TIMES,r,z,mr=200,mz=200)
    one=tight['T']['values'][TABLE_TIMES]
    report['tests']['finite_cylinder_2d_heat']={'max_midplane_vs_1d_K':float(np.max(np.abs(b[:,:,0]-one))),
      'max_end_vs_midplane_K':float(np.max(np.abs(b[:,:,-1]-b[:,:,0]))),
      'series_100_vs_200_max_K':float(np.max(np.abs(a-b))),
      'average_2d_1800':float(mb[-1]),'average_1d_1800':float(tight['T']['average'][-1]),
      'end_center_T_1800':float(b[-1,0,-1]),'end_corner_T_1800':float(b[-1,-1,-1]),
      '2d_mass_executed':False}
    np.savez_compressed(out/'finite_cylinder_2d_heat.npz',time_s=TABLE_TIMES,radius_m=r,z_m=z,T=b,average=mb)
    # Interpolation and transfer coefficient sensitivities; separate from discretization.
    nsens=1280
    lin=cache.get(nsens,{})
    for f in ['T','C']:
        if f not in lin: lin[f]=solve_fv(env,f,n=nsens)
        pp=solve_fv(Environment(data,'pchip'),f,n=nsens)
        report['tests']['pchip_'+f]=differences(lin[f],pp)
        save_field(out/f'pchip_{f}.npz',pp)
    for factor in [.8,1.2]:
        u=solve_fv(env,'C',n=nsens,transfer_scale=factor)
        report['tests'][f'hm_factor_{factor}']={'difference':differences(lin['C'],u),'C1800':u['values'][-1,::5].tolist(),'average1800':float(u['average'][-1])}
        save_field(out/f'hm_{factor}.npz',u)
    # Conditional latent-energy example: intentionally not asserted as calibrated truth.
    le=LatentEnvironment(env,tight['C']['time_s'],tight['C']['values'][:,-1])
    u=solve_fv(le,'T',n=1280,rtol=2e-10,atol=2e-12,max_step=1.,quadrature=True)
    save_field(out/'conditional_latent_T.npz',u)
    V=np.pi*P.R**2*P.L;rho_d=le.rho_d
    water_loss=rho_d*V*(P.C0-tight['C']['average'][-1])
    latent=le.latent*water_loss;sensible=P.rho*P.cp*V*(tight['T']['average'][-1]-P.T0)
    report['tests']['latent_conditional']={'additional_assumptions':['820 interpreted as initial wet bulk density','constant dry density rho/(1+C0)','all effective outward loss evaporates at surface','lambda=2.45 MJ/kg constant (FAO56 Annex3)'],
       'rho_d_kg_dry_m3':rho_d,'lambda_J_kg':le.latent,'implied_water_loss_kg':float(water_loss),
       'baseline_sensible_J':float(sensible),'implied_latent_J':float(latent),'ratio_latent_to_sensible':float(latent/sensible),
       'initial_latent_W_m2':float(le.latent*rho_d*P.hm*(P.C0-data[0,2])),
       'T1800_degC':u['values'][-1,::5].tolist(),'min_T_degC':float(u['values'].min()),
       'max_difference_vs_no_latent_K':float(np.max(np.abs(u['values']-lin['T']['values']))),
       'interpretation':'Diagnostic exposes major closure uncertainty; not used in result1.xlsx; surface trace interpolated at 1s.'}
    report['scales']={'alpha_m2_s':P.alpha,'D_initial_m2_s':diff0,'BiT_R':P.hT*P.R/P.k,
      'BiT_VoverAside':P.hT*P.R/(2*P.k),'BiM_R':P.hm*P.R/diff0,'BiM_VoverAside':P.hm*P.R/(2*diff0),
      'thermal_radial_time_s':P.R**2/P.alpha,'moisture_radial_time_s':P.R**2/diff0,
      'thermal_penetration_m_1800':float(np.sqrt(P.alpha*1800)), 'moisture_penetration_m_1800':float(np.sqrt(diff0*1800)),
      'aspect_ratio_length_diameter':P.L/(2*P.R),'end_to_side_area':P.R/P.L,
      'thermal_center_to_end_time_s':(P.L/2)**2/P.alpha,'moisture_center_to_end_time_s':(P.L/2)**2/diff0}
    report['production_recommendation']={'n':nf,'rtol':2e-12,'atol':2e-14,'max_step_s':3.}
    (out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(report['tests'],indent=2,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
