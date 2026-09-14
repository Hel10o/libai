"""Executable physics-contract, invariance, analytical-reference and sensitivity tests.
Tests are isolated from production inputs/outputs; deliberately faulty operators
are never used to create result2.xlsx.
"""
from __future__ import annotations
import json,time
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import j0,j1,jn_zeros
from q2_core import *
from q2_spectral import Spectral,solve_spectral,spectral_properties

FIELDS=('temperature_degC','moisture_dry_basis')

def difference(a,b,times=None,radii=None):
    d=np.asarray(a)-np.asarray(b);ii=np.unravel_index(np.argmax(abs(d)),d.shape)
    ans={'max_abs':float(np.max(abs(d))),'rms':float(np.sqrt(np.mean(d*d))),'index':list(map(int,ii)),
         'four_decimal_mismatches':int(np.count_nonzero(np.round(a,4)!=np.round(b,4)))}
    if times is not None:ans['time_s']=float(times[ii[0]])
    if radii is not None:ans['radius_cm']=float(radii[ii[1]])
    return ans

def eigen_roots(Bi,n=600):
    hi=jn_zeros(0,n);lo=np.r_[0.,jn_zeros(1,n-1)]
    return np.array([brentq(lambda x:x*j1(x)-Bi*j0(x),a+1e-12,b,xtol=1e-14) for a,b in zip(lo,hi)])

def constant_reference(times,radii,p=P,Te=50.,Ce=.03,nmodes=600):
    S,k,D=properties(np.array([28.]),np.array([2.55]));S,k,D=float(S[0]),float(k[0]),float(D[0])
    ans={}
    for field,diffusivity,Bi,initial,ambient in ((FIELDS[0],k/S,p.hT*p.R/k,p.T0,Te),
                                               (FIELDS[1],D,p.hm*p.R/D,p.C0,Ce)):
        roots=eigen_roots(Bi,nmodes)
        coeff=2*j1(roots)/(roots*(j0(roots)**2+j1(roots)**2))
        shape=j0(roots[:,None]*radii[None,:]/(p.R*100))
        decay=np.exp(-diffusivity*np.asarray(times)[:,None]*roots[None,:]**2/p.R**2)
        vals=ambient+(initial-ambient)*(decay*coeff)@shape
        vals[np.asarray(times)==0]=initial;ans[field]=vals
    return ans

def run_unit_tests(env,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);report={}
    # Dimensional scales from the actual constitutive routines.
    states=[(28,2.55),(50,2.55),(50,.15),(50,.05)]
    report['parameter_states']=[]
    for T,C in states:
        S,k,D=map(float,properties(np.array(T),np.array(C)))
        rho=650+128*C;cp=1450+2736*C/(1+C)
        report['parameter_states'].append({'T_degC':T,'C':C,'rho':rho,'cp':cp,'k':k,'D':D,'S':S,
          'alpha':k/S,'BiT_R':P.hT*P.R/k,'BiC_R':P.hm*P.R/D,'heat_time_R2_alpha_s':P.R**2/(k/S),
          'moisture_time_R2_D_s':P.R**2/D})
    s=report['parameter_states'][0];report['scales']={
      'alpha0':s['alpha'],'D0':s['D'],'BiT_R':s['BiT_R'],'BiC_R':s['BiC_R'],
      'thermal_diffusion_length_3h_m':float(np.sqrt(s['alpha']*10800)),
      'moisture_diffusion_length_3h_m':float(np.sqrt(s['D']*10800)),
      'L_over_diameter':P.L/(2*P.R),'ends_over_side_area':P.R/P.L}
    # Independent 60-digit mpmath evaluation of the printed Appendix-3 formula.
    assert abs(s['D']/5.6416803730256675e-09-1)<2e-14
    # Exact analytical Jacobians checked by complex-step directional differentiation.
    rng=np.random.default_rng(8624)
    op=CoupledFV(40);r=op.r/P.R
    u=np.column_stack([35+3*r*r,2.5-.8*r*r]).ravel();v=rng.normal(size=u.size)
    num=np.imag(op.evaluate(1000,u+1e-28j*v,env))/1e-28
    rel=float(np.max(abs(op.evaluate(1000,u,env,True)@v-num))/np.max(abs(num)))
    sp=Spectral(48);x=sp.x[:-1]
    us=np.column_stack([35+3*x,2.5-.8*x]).ravel();vs=rng.normal(size=us.size)
    nums=np.imag(sp.evaluate(1000,us+1e-28j*vs,env))/1e-28
    rels=float(np.max(abs(sp.evaluate(1000,us,env,True)@vs-nums))/np.max(abs(nums)))
    report['jacobian_complex_step']={'fv_relative':rel,'spectral_relative':rels,'passed':rel<1e-10 and rels<1e-10}
    assert report['jacobian_complex_step']['passed']
    # No driving, actual integration.
    constant=Environment(np.array([[0,28,2.55],[10800,28,2.55]],float))
    z=solve_fv(constant,n=80,end=600,quadrature=True)
    report['uniform_no_drive']={'T_max_change':float(np.max(abs(z[FIELDS[0]]-28))),
                               'C_max_change':float(np.max(abs(z[FIELDS[1]]-2.55)))}
    assert max(report['uniform_no_drive'].values())<1e-10
    # Nonuniform water redistribution with no exchange and exactly uniform temperature.
    closedp=replace(P,hT=0.,hm=0.)
    init=lambda r:np.column_stack([np.full_like(r,35.),2+.4*np.cos(np.pi*(r/P.R)**2)])
    closed=solve_fv(constant,n=160,p=closedp,end=600,initial=init,quadrature=True)
    report['closed_redistribution']={'T_max_change':float(np.max(abs(closed[FIELDS[0]]-35))),
      'C_weighted_average_drift':float(np.max(abs(closed['average_moisture_dry_basis']-closed['average_moisture_dry_basis'][0]))),
      'water_profile_changed':float(np.max(abs(closed[FIELDS[1]][-1]-closed[FIELDS[1]][0]))),
      'scope':'Fixed-reference-dry-density effective mass integral; not rho(C)*C'}
    assert report['closed_redistribution']['T_max_change']<1e-8
    assert report['closed_redistribution']['C_weighted_average_drift']<1e-9
    # Frozen coefficients versus cylindrical Bessel solution.
    bench=Environment(np.array([[0,50,.03],[600,50,.03]],float))
    analytical_sol=solve_spectral(bench,n=128,end=600,freeze=True,rtol=1e-12)
    analytic=constant_reference(analytical_sol['time_s'],analytical_sol['radius_cm'])
    report['constant_coefficient_bessel']={k:difference(analytical_sol[k][1:],analytic[k][1:],analytical_sol['time_s'][1:],analytical_sol['radius_cm']) for k in FIELDS}
    assert report['constant_coefficient_bessel'][FIELDS[0]]['max_abs']<2e-7
    assert report['constant_coefficient_bessel'][FIELDS[1]]['max_abs']<2e-8
    save_solution(out/'constant_coefficient_spectral.npz',analytical_sol)
    np.savez_compressed(out/'constant_coefficient_bessel.npz',**analytic,time_s=analytical_sol['time_s'])
    # Fault injection: these tests must identify incorrect model implementations.
    faults={}
    correctD=s['D'];wrongD=.0024*np.exp(-.45/2.55-3850/28)
    faults['Celsius_in_exponent']={'wrong_over_correct_D':float(wrongD/correctD),'detected':abs(wrongD/correctD-1)>.99}
    mid=(r>.1)&(r<.9);expected=4/P.R**2/properties(u[0::2],u[1::2])[0]*3*(
      properties(u[0::2],u[1::2])[1]+r*r*.38/(1+u[1::2])**2*(-.8))
    correct=op.evaluate(1000,u,env)[0::2]
    faulty=CoupledFV(40,fault='face_alpha').evaluate(1000,u,env)[0::2]
    ec=float(np.max(abs(correct[mid]-expected[mid])));ef=float(np.max(abs(faulty[mid]-expected[mid])))
    faults['face_alpha_instead_of_k']={'correct_interior_error_K_s':ec,'faulty_error_K_s':ef,'detected':ef>100*ec}
    heated=u.copy();heated[0::2]+=10
    reaction=np.max(abs(op.evaluate(1000,heated,env)[1::2]-op.evaluate(1000,u,env)[1::2]))
    fr=CoupledFV(40,fault='frozen_D')
    reaction_bad=np.max(abs(fr.evaluate(1000,heated,env)[1::2]-fr.evaluate(1000,u,env)[1::2]))
    faults['frozen_RHS_D_coupling']={'correct_response_C_s':float(reaction),'faulty_response':float(reaction_bad),'detected':reaction>1e-6 and reaction_bad<1e-20}
    uniform=init(op.r).ravel()
    true_rhs=CoupledFV(40,p=closedp).evaluate(0,uniform,constant)[0::2]
    bad_rhs=CoupledFV(40,p=closedp,fault='product_storage').evaluate(0,uniform,constant)[0::2]
    faults['unclosed_product_storage']={'correct_T_rate_max':float(np.max(abs(true_rhs))),
       'faulty_T_rate_max':float(np.max(abs(bad_rhs))),'detected':np.max(abs(true_rhs))<1e-14 and np.max(abs(bad_rhs))>1e-5}
    report['fault_injection']=faults
    assert all(v['detected'] for v in faults.values())
    # Unit change, separate integrations on the same mesh.
    celsius=solve_fv(env,n=320,rtol=3e-13,atol_T=2e-14,atol_C=2e-15,max_step=10.)
    kelvin=solve_fv(env,n=320,rtol=3e-13,atol_T=2e-14,atol_C=2e-15,max_step=10.,temperature_offset=273.15)
    report['temperature_unit_invariance']={k:difference(celsius[k],kelvin[k],celsius['time_s'],celsius['radius_cm']) for k in FIELDS}
    assert report['temperature_unit_invariance'][FIELDS[0]]['max_abs']<2e-8
    assert report['temperature_unit_invariance'][FIELDS[1]]['max_abs']<2e-9
    save_solution(out/'kelvin_invariance.npz',kelvin)
    report['passed']=True
    (out/'unit_tests.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,default=lambda x:x.item()),encoding='utf-8')
    return json.loads(json.dumps(report,default=lambda x:x.item()))

def run_sensitivity(env,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    baseline=solve_fv(env,n=320,rtol=1e-10,max_step=20.)
    save_solution(out/'sensitivity_baseline.npz',baseline);report={}
    cases=[('hT_0.8',replace(P,hT=.8*P.hT),env,False),('hT_1.2',replace(P,hT=1.2*P.hT),env,False),
           ('hm_0.8',replace(P,hm=.8*P.hm),env,False),('hm_1.2',replace(P,hm=1.2*P.hm),env,False),
           ('pchip_all241',P,Environment(env.data,kind='pchip'),False),('freeze_all_initial',P,env,True)]
    for name,p,e,freeze in cases:
        sol=solve_fv(e,n=320,p=p,freeze=freeze,rtol=1e-10,max_step=20.)
        report[name]={'differences':{k:difference(sol[k],baseline[k],sol['time_s'],sol['radius_cm']) for k in FIELDS},
          'T_center_end':float(sol[FIELDS[0]][-1,0]),'T_surface_end':float(sol[FIELDS[0]][-1,-1]),
          'C_center_end':float(sol[FIELDS[1]][-1,0]),'C_surface_end':float(sol[FIELDS[1]][-1,-1]),
          'mean_C_end':float(sol['average_moisture_dry_basis'][-1]),'stats':sol['stats']}
        save_solution(out/f'sensitivity_{name}.npz',sol)
    report['scope']='Deterministic scenarios, not confidence intervals; all compared with identical-mesh baseline'
    (out/'sensitivity.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,default=lambda x:x.item()),encoding='utf-8')
    return json.loads(json.dumps(report,default=lambda x:x.item()))
