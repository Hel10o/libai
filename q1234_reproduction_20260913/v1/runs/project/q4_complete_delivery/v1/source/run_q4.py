"""Cold-start fully coupled Q4 material-coordinate calculation.

The effective equations are S(C) D_s T=div(k grad T), D_s C=div(D grad C).
With v_s=w=(R'/R)r, x=(r/R)^2, radial div=4/R² d_x(x flux_x).
No mass-density dilution term is added to dry-basis C. Rho(C) is only in S.
"""
import argparse,json,time,hashlib,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from numpy.polynomial.chebyshev import chebvander
from q4_common import ROOT,Environment,Radius,read_numeric,load_input,input_audit
from q4_spectral import Reference,interp,primitive,inverse_primitive

def run(n,out,kind='linear',method='Radau',rtol=2e-10,max_step=180.,audit=False):
    out=Path(out)
    if out.with_suffix('.json').exists() or out.with_suffix('.npz').exists():raise FileExistsError(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    data=load_input(ROOT/'inputs/attachment1.xlsx');rd=read_numeric(ROOT/'inputs/attachment2.xlsx',2)[1]
    env=Environment(data);radius=Radius(rd,kind);op=Reference(n,radius)
    weights=np.linalg.solve(chebvander(2*op.x-1,n).T,np.array([0. if j%2 else 1/(1-j*j) for j in range(n+1)]))
    target=np.linspace(0,.02,21);times=[0.];states=[np.tile([28.,2.55],(n+1,1))]
    radii=[.02];fixed=[np.tile([28.,2.55],(21,1))];surf=[[28.,2.55]];meanC=[2.55]
    maxima=[2.55];positions=[0.];rate_checks=[];cumulative_out=0.;balance=[]
    y=np.tile([28.,2.55],n);tic=time.perf_counter();steps=nfev=nlu=0;event_info=None
    cuts=np.unique(np.r_[data[:,0],rd[:,0],np.arange(259200,518401,1800.)]);next_out=60.
    gauss_x,gauss_w=np.polynomial.legendre.leggauss(4)
    def record(t,u,e):
        R=radius(t);inside=target<=R+1e-14;v=np.full((21,2),np.nan)
        I=interp(op.x,op.w,(target[inside]/R)**2)
        v[inside,0]=I@u[:,0];v[inside,1]=inverse_primitive(I@primitive(u[:,1]),I@u[:,1])
        ex=op.extrema(u)
        times.append(float(t));states.append(u);radii.append(R);fixed.append(v);surf.append(u[-1]);meanC.append(float(weights@u[:,1]));maxima.append(ex['max_C']);positions.append(float(R*np.sqrt(ex['x_max'])))
    for a,b in zip(cuts[:-1],cuts[1:]):
        e=env.segment(a,b);fun=lambda t,y:op.evaluate(t,y,e)
        def ev(t,y):return float(np.max(y[1::2])-.15)
        ev.direction=-1;ev.terminal=False
        sol=solve_ivp(fun,(a,b),y,method=method,jac=lambda t,y:op.evaluate(t,y,e,True),
             rtol=rtol,atol=np.tile([rtol*.1,rtol*.01],n),max_step=min(max_step,30.) if a<14400 else max_step,dense_output=True,events=ev)
        if not sol.success:raise RuntimeError(sol.message)
        steps+=len(sol.t)-1;nfev+=sol.nfev;nlu+=sol.nlu;end=b
        if len(sol.t_events[0]):
            root=float(sol.t_events[0][0]);end=.36*np.ceil(root/.36)
            # Exact evaluations on an explicitly stored grid allow the final
            # execution margin to be selected AFTER comparing numerical runs.
            near=np.unique(np.r_[root+np.array([-2.,-1.,0.,1.,2.]),
                      np.arange(np.floor(root/.36)-3,np.ceil(root/.36)+10)*.36])
            if near.min()<a or near.max()>b:raise RuntimeError('Event neighborhood crosses segment; refine cuts')
            near_full=np.array([op.surface(t,sol.sol(t).reshape(n,2),e) for t in near])
            event_info={'critical_s':root,'critical_h':root/3600,'preliminary_execution_s':end,
                 'preliminary_execution_max_C':op.extrema(op.surface(end,sol.sol(end).reshape(n,2),e))['max_C'],
                 'critical_slope':float(fun(root,sol.sol(root))[1]),'root_extrema':op.extrema(op.surface(root,sol.sol(root).reshape(n,2),e)),
                 'radius_m':radius(root),'length_m':.25,'root_is_1d':True}
        ts=np.arange(next_out,end+1e-9,60.)
        if event_info:ts=np.r_[ts[ts<end-1e-8],end]
        for t in ts:record(t,op.surface(t,sol.sol(t).reshape(n,2),e),e)
        if len(ts):next_out=60*(np.floor(ts[-1]/60)+1)
        # Independent quadrature of the boundary loss; dry mass cancels.
        if audit:
            for sa,sb in zip(sol.t[:-1],sol.t[1:]):
                sb=min(sb,end)
                if sa>=end:break
                tq=(sa+sb)/2+(sb-sa)/2*gauss_x
                rates=[]
                for t in tq:
                    u=op.surface(t,sol.sol(t).reshape(n,2),e)
                    rates.append(2/radius(t)*8e-7*(u[-1,1]-e(t)[1]))
                cumulative_out+=(sb-sa)/2*np.dot(gauss_w,rates)
            u=op.surface(end,sol.sol(end).reshape(n,2),e)
            balance.append([float(end),float(weights@u[:,1]),float(cumulative_out),float(weights@u[:,1]+cumulative_out-2.55)])
        if b%21600==0 or event_info:
            u=op.surface(end,sol.sol(end).reshape(n,2),e)
            f=fun(end,sol.sol(end)).reshape(n,2)
            # Pointwise differentiated polynomial at the boundary supplies rate
            # solely for this spatial quadrature identity, not a PDE unknown.
            T,C=u.T;R=radius(end);S=(760+90*C)*(1850+2150*C/(1+C));k=.12+.20*C/(1+C)
            U=primitive(C);gt=op.G@(T-T[-1]);gu=op.G@(U-U[-1]);qt=k*gt;qm=.00042*np.exp(-3850/(T+273.15))*gu
            ft=4/R**2*(qt+op.x*(op.G@(qt-qt[-1])))/S;fc=4/R**2*(qm+op.x*(op.G@(qm-qm[-1])))
            rate_checks.append({'t_s':end,'water_rate_residual':float(weights@fc+2/R*8e-7*(C[-1]-e(end)[1])),
                'effective_heat_rate_residual':float(weights@(S*ft)+2/R*25*(T[-1]-e(end)[0]))})
            print(json.dumps({'n':n,'kind':kind,'time_h':end/3600,'max_C':float(C.max()),'elapsed_s':time.perf_counter()-tic}),flush=True)
        y=sol.sol(end)
        if event_info:break
    if event_info is None:raise RuntimeError('No threshold by 144h; no fabricated endpoint')
    stats={'n':n,'method':method,'rtol':rtol,'atol_TC':[rtol*.1,rtol*.01],'max_step_s':max_step,'radius_interpolation':kind,
       'future_environment':env.future.tolist(),'event':event_info,'elapsed_s':time.perf_counter()-tic,'steps':steps,'nfev':nfev,'nlu':nlu,
       'max_boundary_residual':op.max_bc_residual,'boundary_fallbacks':op.boundary_fallbacks,'rate_checks':rate_checks,
       'mass_balance_max_abs_in_mean_C':float(np.max(np.abs(np.array(balance)[:,-1]))) if balance else None,
       'numpy':np.__version__,'scipy':scipy.__version__,'python':platform.python_version(),
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    np.savez_compressed(out.with_suffix('.npz'),time_s=times,radius_m=radii,full_TC=states,x=op.x,
          fixed_radius_m=target,sample_TC=fixed,surface_TC=surf,mean_C=meanC,max_C=maxima,max_radius_m=positions,
          event_time_s=near,event_full_TC=near_full,quadrature_weights=weights,balance=balance)
    out.with_suffix('.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2))
    print(json.dumps(stats,ensure_ascii=False),flush=True)
    return stats

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=80);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--kind',choices=['linear','fixed','pchip'],default='linear');p.add_argument('--method',default='Radau')
    p.add_argument('--rtol',type=float,default=2e-10);p.add_argument('--max-step',type=float,default=180.);p.add_argument('--audit',action='store_true')
    a=p.parse_args();run(a.n,a.out,a.kind,a.method,a.rtol,a.max_step,a.audit)
