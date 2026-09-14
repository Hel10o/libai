"""Cold start a unified, simultaneous Appendix 3 T,C trajectory at every second.

The copied historical Chebyshev primitive-flux spatial implementation is used;
its floating-point output is retained, never cell-patched or spliced. Each data
knot is an integration boundary. After 4 h the same Q3 arithmetic mean applies.
Vectorization below changes only dense-output profile evaluation, not the ODE.
"""
from pathlib import Path
import sys, json, time, argparse, hashlib, platform
sys.path.insert(0,str(Path(__file__).resolve().parent/'legacy_q3'))
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from q3_reference import Reference, primitive, inverse_primitive
from q3_solver import Environment, load_input

def batched_profiles(op, ts, ys, env):
    # ys = (samples,n,2); solve one Robin surface per sample.
    Ti,Ci=ys[:,:,0],ys[:,:,1]
    et,ec=np.asarray(env(ts)).T
    Ui=primitive(Ci);anchor=Ci[:,-1];ua=Ui[:,-1]
    gt=(Ti-et[:,None])@op.row;gu=(Ui-ua[:,None])@op.row
    c=anchor.copy()
    for it in range(40):
        k=.21+.38*c/(1+c);kc=.38/(1+c)**2
        den=op.bd*k*op.diag+25
        Ts=et-op.bd*k*gt/den
        Tsc=-op.bd*kc*gt*25/den**2
        A=.0024*np.exp(-3850/(Ts+273.15));At=A*3850/(Ts+273.15)**2
        ux=gu+op.diag*(primitive(c)-ua)
        val=op.bd*A*ux+8e-7*(c-ec)
        deriv=op.bd*(At*Tsc*ux+A*op.diag*np.exp(-.45/c))+8e-7
        step=val/deriv
        bad=c-step<=0
        while bad.any(): step[bad]*=.5;bad=c-step<=0
        c-=step
        if np.max(abs(step)/(1+abs(c)))<2e-14:break
    else: raise RuntimeError('Vectorized Robin Newton failed')
    k=.21+.38*c/(1+c);Ts=et-op.bd*k*gt/(op.bd*k*op.diag+25)
    full=np.concatenate([ys,np.stack([Ts,c],axis=1)[:,None,:]],axis=1)
    T,C=full[:,:,0],full[:,:,1]
    U=primitive(C)
    profile=np.stack([T@op.I.T,inverse_primitive(U@op.I.T,C@op.I.T)],axis=-1)
    return profile,full

def solve(n=160,end=206906.76,method='Radau',rtol=2e-12,max_step=60.,out=None,input_path=None,step=1.):
    out=Path(out)
    if out.with_suffix('.npz').exists() or out.with_suffix('.json').exists():raise FileExistsError(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    op=Reference(n);data=load_input(input_path);env=Environment(data,'mean','linear')
    ts=np.unique(np.r_[np.arange(0,end+1e-8,step),end]);v=np.empty((len(ts),21,2));v[0]=[28.,2.55]
    y=np.tile([28.,2.55],n)
    cuts=np.unique(np.r_[data[:,0][data[:,0]<end],np.arange(21600.,end,21600.),end]);snapt=[0.];snaps=[np.tile([28.,2.55],(n+1,1))]
    stats=dict(n=n,method=method,rtol=rtol,atol_T=2e-12,atol_C=2e-14,max_step_s=max_step,end_s=end,step_s=step,source_sha='anonymous-source',model='Appendix 3 effective radial, no explicit latent heat; same equation and boundary as Q2/Q3',environment='piecewise linear through 4 h; then arithmetic mean of 61 samples from 3 to 4 h',future=env.future.tolist(),nsteps=0,nfev=0,nlu=0,events=[],runtime_versions={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__})
    start=time.perf_counter();max_batch_diff=np.zeros(2)
    for a,b in zip(cuts[:-1],cuts[1:]):
        ee=env.segment(a,b)
        def event(t,y):return y[1::2].max()-.15
        event.terminal=False;event.direction=-1
        sol=solve_ivp(lambda t,y:op.evaluate(t,y,ee),(a,b),y,method=method,jac=lambda t,y:op.evaluate(t,y,ee,True),rtol=rtol,atol=np.tile([2e-12,2e-14],n),max_step=min(max_step,20.) if a<14400 else max_step,dense_output=True,events=event)
        if not sol.success:raise RuntimeError(sol.message)
        stats['nsteps']+=len(sol.t)-1;stats['nfev']+=sol.nfev;stats['nlu']+=sol.nlu
        ids=np.flatnonzero((ts>a)&(ts<=b+1e-8))
        for chunk in np.array_split(ids,max(1,(len(ids)+2047)//2048)):
            if len(chunk)==0:continue
            yy=sol.sol(ts[chunk]).reshape(n,2,-1).transpose(2,0,1)
            vv,full=batched_profiles(op,ts[chunk],yy,ee);v[chunk]=vv
            # Scalar evaluation samples ensure vectorization has not changed equations.
            for j in [0,len(chunk)-1]:
                scalar=op.profiles(op.surface(ts[chunk[j]],yy[j],ee))
                max_batch_diff=np.maximum(max_batch_diff,np.max(abs(scalar-vv[j]),axis=0))
        for t in sol.t_events[0]:
            state=op.surface(t,sol.sol(t).reshape(n,2),ee)
            stats['events'].append({'root_s':float(t),'max_C':float(state[:,1].max()),'extrema':op.extrema(state)})
        y=sol.y[:,-1]
        if b%1800==0 or b==end:
            sf=op.surface(b,y.reshape(n,2),ee);snapt.append(b);snaps.append(sf)
            print(json.dumps({'n':n,'t_s':float(b),'Cmax':float(sf[:,1].max()),'elapsed_s':time.perf_counter()-start}),flush=True)
    stats['elapsed_s']=time.perf_counter()-start;stats['vectorized_vs_scalar_max_abs_TC']=max_batch_diff.tolist();stats['end_max_C']=float(snaps[-1][:,1].max());stats['end_global_extrema']=op.extrema(snaps[-1]);stats['strict_1d_pass']=stats['end_global_extrema']['max_C']<.15
    stats['input_sha256']=hashlib.sha256(Path(input_path).read_bytes()).hexdigest();stats['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    np.savez_compressed(out.with_suffix('.npz'),time_s=ts,radius_cm=np.linspace(0,2,21),profile_TC=v,environment=env(ts),x=op.x,snapshot_time_s=np.array(snapt),snapshot_full_TC=np.array(snaps))
    out.with_suffix('.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(stats,ensure_ascii=False),flush=True)
    return stats

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=160);p.add_argument('--end',type=float,default=206906.76);p.add_argument('--method',default='Radau');p.add_argument('--rtol',type=float,default=2e-12);p.add_argument('--max-step',type=float,default=60.);p.add_argument('--step',type=float,default=1.);p.add_argument('--out',type=Path,required=True);p.add_argument('--input',type=Path,default=Path(__file__).resolve().parents[1]/'inputs/attachment1.xlsx');a=p.parse_args();solve(a.n,a.end,a.method,a.rtol,a.max_step,a.out,a.input,a.step)
