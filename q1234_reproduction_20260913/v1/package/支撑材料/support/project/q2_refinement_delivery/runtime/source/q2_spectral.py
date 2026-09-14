"""Independent Chebyshev collocation in x=(r/R)^2.
No FV operator or face averaging is reused. Robin conditions are algebraically
eliminated. All nonlinear material derivatives and boundary derivatives are
included in the dense Jacobian. BDF/Radau integrate the interior system.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from q2_core import P,Environment,load_environment,save_solution

def spectral_properties(T,C,freeze=False):
    if np.any(np.real(C)<=0) or np.any(np.real(T+273.15)<=0):raise FloatingPointError('Invalid spectral state')
    if freeze:T=T*0+28.;C=C*0+2.55
    rho=650.+128.*C
    cp=(1450.+4186.*C)/(1.+C)
    k=(.21+.59*C)/(1.+C)
    S=rho*cp;D=.0024*np.exp(-.45/C)*np.exp(-3850./(T+273.15))
    Sc=128.*cp+rho*2736./(1.+C)**2
    kc=.38/(1.+C)**2;Dt=3850.*D/(T+273.15)**2;Dc=.45*D/C**2
    if freeze:Sc*=0;kc*=0;Dt*=0;Dc*=0
    return S,k,D,Sc,kc,Dt,Dc

def cheb_grid(n):
    j=np.arange(n+1);x=(1-np.cos(np.pi*j/n))/2
    bw=(-1.)**j;bw[[0,-1]]*=.5
    dx=x[:,None]-x[None,:]
    G=np.divide(bw[None,:]/bw[:,None],dx,out=np.zeros_like(dx),where=dx!=0)
    G[np.diag_indices(n+1)]=-G.sum(axis=1)
    # Derivative implementation anchors at the surface to preserve constants.
    G[:,-1]=-G[:,:-1].sum(axis=1)
    theta=np.pi*j/n;weights=np.zeros(n+1);ii=np.arange(1,n);v=np.ones(n-1)
    if n%2==0:
        weights[[0,n]]=1/(n*n-1)
        for k in range(1,n//2):v-=2*np.cos(2*k*theta[ii])/(4*k*k-1)
        v-=np.cos(n*theta[ii])/(n*n-1)
    else:
        weights[[0,n]]=1/(n*n)
        for k in range(1,(n-1)//2+1):v-=2*np.cos(2*k*theta[ii])/(4*k*k-1)
    weights[ii]=2*v/n
    return x,G,bw,weights/2

def interpolation_matrix(x,bw,target):
    M=np.empty((len(target),len(x)))
    for i,t in enumerate(target):
        d=t-x;j=np.argmin(abs(d))
        if abs(d[j])<1e-14:M[i]=0.;M[i,j]=1.
        else:
            v=bw/d;M[i]=v/v.sum()
    return M

class Spectral:
    def __init__(self,n,p=P,freeze=False):
        self.n=n;self.m=n+1;self.p=p;self.freeze=freeze
        self.x,self.G,self.bw,self.weights=cheb_grid(n)
        self.row=self.G[-1,:-1];self.diag=-self.row.sum();self.fac=4/p.R**2;self.A=2/p.R
        self.interp=interpolation_matrix(self.x,self.bw,np.linspace(0,1,21)**2)
        self.max_newton=0
    def surface(self,t,u,env,derivative=False):
        Ti=u[:,0];Ci=u[:,1];et,ec=env(t)
        gT=self.row@(Ti-et);anchor=Ci[-1];gC=self.row@(Ci-anchor)
        near=anchor-gC/self.diag
        def calc(c,need=False):
            _,k,_,_,kc,_,_=spectral_properties(et+0*c,c,self.freeze)
            den=self.A*k*self.diag+self.p.hT
            T=et-self.A*k*gT/den
            Tp=-self.A*kc*gT*self.p.hT/den**2
            S,k,D,Sc,kc,Dt,Dc=spectral_properties(T,c,self.freeze)
            cx=gC+self.diag*(c-anchor)
            F=self.A*D*cx+self.p.hm*(c-ec)
            Fp=self.A*((Dc+Dt*Tp)*cx+D*self.diag)+self.p.hm
            return F,Fp,T,cx,(S,k,D,Sc,kc,Dt,Dc)
        c=near
        for it in range(20):
            F,Fp,T,cx,props=calc(c)
            step=F/Fp;trial=c-step
            while np.real(trial)<=0:step*=.5;trial=c-step
            c=trial
            if abs(step)<3e-14*(1+abs(c)):break
        else:
            # On an unexpected Newton failure choose the largest positive root
            # continuously connected to the Neumann-side value, never a tiny spurious root.
            if np.iscomplexobj(u):raise RuntimeError('Complex-step surface solve failed')
            high=float(near);fh=calc(high)[0];found=False
            for low in np.geomspace(high*.99,max(1e-12,high*1e-9),160):
                fl=calc(low)[0]
                if fl*fh<=0:
                    c=brentq(lambda a:calc(a)[0],low,high,xtol=5e-15);found=True;break
                high=low;fh=fl
            if not found:raise RuntimeError('Positive surface root not found')
        self.max_newton=max(self.max_newton,it+1)
        F,Fp,T,cx,props=calc(c)
        full=np.vstack([u,np.array([T,c])])
        if not derivative:return full
        _,k,D,_,kc,Dt,Dc=props
        tx=self.row@(Ti-T)
        M=np.array([[self.A*k*self.diag+self.p.hT,self.A*kc*tx],
                    [self.A*Dt*cx,self.A*(D*self.diag+Dc*cx)+self.p.hm]])
        part=np.zeros((2,2*self.n));part[0,0::2]=self.A*k*self.row;part[1,1::2]=self.A*D*self.row
        bd=-np.linalg.solve(M,part)
        BT=np.zeros((self.m,2*self.n));BC=np.zeros_like(BT)
        BT[np.arange(self.n),2*np.arange(self.n)]=1.
        BC[np.arange(self.n),2*np.arange(self.n)+1]=1.
        BT[-1]=bd[0];BC[-1]=bd[1]
        return full,BT,BC
    def evaluate(self,t,y,env,jac=False):
        u=y.reshape(self.n,2)
        if jac:full,BT,BC=self.surface(t,u,env,True)
        else:full=self.surface(t,u,env)
        T,C=full.T;S,k,D,Sc,kc,Dt,Dc=spectral_properties(T,C,self.freeze)
        gt=self.G@(T-T[-1]);gc=self.G@(C-C[-1])
        qt=k*gt;qc=D*gc
        ht=qt+self.x*(self.G@(qt-qt[-1]));hc=qc+self.x*(self.G@(qc-qc[-1]))
        fT=self.fac*ht/S;fC=self.fac*hc
        if not jac:return np.column_stack([fT[:-1],fC[:-1]]).ravel()
        gtj=self.G@(BT-BT[-1]);gcj=self.G@(BC-BC[-1])
        tqj=k[:,None]*gtj+(kc*gt)[:,None]*BC
        cqj=D[:,None]*gcj+gc[:,None]*(Dt[:,None]*BT+Dc[:,None]*BC)
        htj=tqj+self.x[:,None]*(self.G@(tqj-tqj[-1]))
        hcj=cqj+self.x[:,None]*(self.G@(cqj-cqj[-1]))
        ftj=self.fac*htj/S[:,None]-(fT*Sc/S)[:,None]*BC
        fcj=self.fac*hcj
        J=np.empty((2*self.n,2*self.n));J[0::2]=ftj[:-1];J[1::2]=fcj[:-1]
        return J

def solve_spectral(env,n=160,p=P,end=10800.,rtol=2e-12,atol_T=2e-13,atol_C=2e-14,max_step=15.,method='BDF',freeze=False,verbose=False):
    op=Spectral(n,p,freeze);y=np.tile([p.T0,p.C0],n)
    times=np.arange(int(end)+1,dtype=float)
    if times[-1]!=end:raise ValueError('Spectral test end must be integer seconds')
    out=np.empty((len(times),21,2));out[0,:,0]=p.T0;out[0,:,1]=p.C0
    average=np.empty((len(times),2));average[0]=[p.T0,p.C0]
    snapshot_times=[0.];snaps=[np.tile([p.T0,p.C0],(n+1,1))]
    cuts=np.unique(np.r_[0.,env.knots[(env.knots>0)&(env.knots<end)],end])
    stats={'n':n,'spatial_method':'Chebyshev-Lobatto collocation in x=(r/R)^2; coupled Robin elimination',
      'time_method':method,'rtol':rtol,'atol_T':atol_T,'atol_C':atol_C,'max_step_s':max_step,
      'nfev':0,'njev':0,'nlu':0,'accepted_steps':0,'environment':env.kind,'extension':env.extension,'end_s':end}
    tic=time.perf_counter()
    for a,b in zip(cuts[:-1],cuts[1:]):
        sol=solve_ivp(lambda t,y:op.evaluate(t,y,env),(a,b),y,method=method,
          jac=lambda t,y:op.evaluate(t,y,env,True),rtol=rtol,atol=np.tile([atol_T,atol_C],n),
          max_step=max_step,dense_output=True,first_step=min(1e-4,b-a))
        if not sol.success:raise RuntimeError(f'Spectral failed on [{a},{b}]: {sol.message}')
        for k in ('nfev','njev','nlu'):stats[k]+=int(getattr(sol,k))
        stats['accepted_steps']+=len(sol.t)-1
        mask=(times>a)&(times<=b+1e-8);ts=times[mask];ys=sol.sol(ts).reshape(n,2,-1).transpose(2,0,1)
        for ix,t,u in zip(np.flatnonzero(mask),ts,ys):
            full=op.surface(t,u,env);out[ix]=op.interp@full;average[ix]=op.weights@full
            if t in (1800.,3600.,5400.,7200.,9000.,10800.):snapshot_times.append(t);snaps.append(full.copy())
        if not np.isfinite(out[mask]).all() or np.any(out[mask,:,1]<=0):raise FloatingPointError('Spectral output invalid')
        y=sol.y[:,-1]
        if verbose and int(b)%1800==0:print(f'Cheb n={n}, t={b:g}',flush=True)
    stats['elapsed_s']=time.perf_counter()-tic;stats['max_boundary_newton_iterations']=op.max_newton
    return {'time_s':times,'radius_cm':np.linspace(0,2,21),'temperature_degC':out[:,:,0],
      'moisture_dry_basis':out[:,:,1],'average_temperature_degC':average[:,0],
      'average_moisture_dry_basis':average[:,1],'environment':env(times),
      'internal_radius_m':p.R*np.sqrt(op.x),'internal_weights_normalized':op.weights,
      'snapshot_times_s':np.array(snapshot_times),'snapshots_temperature_degC':np.array(snaps)[:,:,0],
      'snapshots_moisture_dry_basis':np.array(snaps)[:,:,1],'stats':stats}

def main():
    b=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=b/'inputs/附件1.xlsx');ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n',type=int,default=160);ap.add_argument('--rtol',type=float,default=2e-12)
    ap.add_argument('--end',type=float,default=10800.);ap.add_argument('--method',choices=['BDF','Radau'],default='BDF')
    args=ap.parse_args();a,_=load_environment(args.input)
    sol=solve_spectral(Environment(a),n=args.n,end=args.end,rtol=args.rtol,method=args.method,verbose=True)
    save_solution(args.out/f'cheb_n{args.n}_{args.method}.npz',sol);print(json.dumps(sol['stats'],indent=2))
if __name__=='__main__':main()
