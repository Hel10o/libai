"""Independent global Chebyshev collocation with primitive-gradient moisture flux.
No FV control volumes, face weights, or FV spatial operator are reused.
Full coupled heat/moisture evolution from t=0; Robin boundary algebraically solved.
Grid/differentiation construction follows the audited q2_spectral.py approach;
moisture flux and its boundary Jacobian are newly derived for long-time Q3.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from scipy.special import exp1
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.fft import dct
from q3_solver import load_input,Environment

def primitive(c):
 if np.any(np.real(c)<=0):raise FloatingPointError('C<=0 in reference')
 return c*np.exp(-.45/c)-.45*exp1(.45/c)

def grid(n):
 x=(1-np.cos(np.pi*np.arange(n+1)/n))/2
 w=(-1.)**np.arange(n+1);w[[0,-1]]*=.5
 dx=x[:,None]-x[None,:]
 G=np.divide(w[None,:]/w[:,None],dx,out=np.zeros_like(dx),where=dx!=0)
 G[np.diag_indices(n+1)]=-G.sum(axis=1);G[:,-1]=-G[:,:-1].sum(axis=1)
 return x,G,w

def interp(x,w,target):
 d=target[:,None]-x
 out=np.empty_like(d)
 for j,dd in enumerate(d):
  k=np.argmin(abs(dd))
  if abs(dd[k])<1e-14:out[j]=0;out[j,k]=1
  else:v=w/dd;out[j]=v/v.sum()
 return out

def inverse_primitive(v,guess):
 c=np.asarray(guess).copy()
 for _ in range(40):
  delta=(primitive(c)-v)/np.exp(-.45/c)
  cnew=c-delta
  # Globalize Newton without altering the equation or accepted value.
  bad=np.real(cnew)<=0
  while np.any(bad):delta[bad]*=.5;cnew=c-delta;bad=np.real(cnew)<=0
  c=cnew
  if np.max(abs(delta))<3e-14:break
 else:raise RuntimeError('Inverse primitive failed')
 return c

class Reference:
 def __init__(self,n):
  self.n=n;self.x,self.G,self.w=grid(n);self.row=self.G[-1,:-1];self.diag=self.G[-1,-1]
  self.fac=4/.02**2;self.bd=2/.02
  self.BT0=np.zeros((n+1,2*n));self.BC0=np.zeros_like(self.BT0)
  self.BT0[np.arange(n),2*np.arange(n)]=1;self.BC0[np.arange(n),2*np.arange(n)+1]=1
  self.I=interp(self.x,self.w,np.linspace(0,1,21)**2)
  self.denseI=interp(self.x,self.w,np.linspace(0,1,8*n+1))
  self.min_boundary_derivative=float('inf');self.max_bc_residual=0.;self.max_newton=0;self.boundary_fallbacks=0
 def surface(self,t,u,env,jac=False):
  Ti,Ci=u.T;et,ec=env(t);Ui=primitive(Ci);anchor=Ci[-1];ua=Ui[-1]
  gt=self.row@(Ti-et);gu=self.row@(Ui-ua)
  def calc(c):
   k=.21+.38*c/(1+c);kc=.38/(1+c)**2
   den=self.bd*k*self.diag+25
   Ts=et-self.bd*k*gt/den
   Tsc=-self.bd*kc*gt*25/den**2
   A=.0024*np.exp(-3850/(Ts+273.15));At=A*3850/(Ts+273.15)**2
   ux=gu+self.diag*(primitive(c)-ua)
   f=np.exp(-.45/c)
   val=self.bd*A*ux+8e-7*(c-ec)
   deriv=self.bd*(At*Tsc*ux+A*self.diag*f)+8e-7
   return val,deriv,Ts,k,kc,A,At,ux
  c=anchor
  for it in range(25):
   val,deriv,*_=calc(c);step=val/deriv
   while c-step<=0:step*=.5
   c=c-step
   if abs(step)<2e-14*(1+abs(c)):break
  else:
   self.boundary_fallbacks+=1
   lo=min(ec,float(Ci.min()))*.1;hi=max(ec,float(Ci.max()))*1.1
   c=brentq(lambda z:calc(z)[0],lo,hi,xtol=5e-15)
  val,deriv,Ts,k,kc,A,At,ux=calc(c)
  self.max_newton=max(self.max_newton,it+1);self.min_boundary_derivative=min(self.min_boundary_derivative,float(deriv));self.max_bc_residual=max(self.max_bc_residual,float(abs(val)))
  full=np.vstack([u,[Ts,c]])
  if not jac:return full
  tx=self.row@(Ti-Ts)
  B=np.array([[self.bd*k*self.diag+25,self.bd*kc*tx],[self.bd*At*ux,self.bd*A*self.diag*np.exp(-.45/c)+8e-7]])
  Q=np.zeros((2,2*self.n));Q[0,::2]=self.bd*k*self.row;Q[1,1::2]=self.bd*A*self.row*np.exp(-.45/Ci)
  sens=-np.linalg.solve(B,Q);BT=self.BT0.copy();BC=self.BC0.copy();BT[-1]=sens[0];BC[-1]=sens[1]
  return full,BT,BC
 def evaluate(self,t,y,env,jac=False):
  if jac:full,BT,BC=self.surface(t,y.reshape(self.n,2),env,True)
  else:full=self.surface(t,y.reshape(self.n,2),env)
  T,C=full.T
  if np.any(T<=-273.15):raise FloatingPointError('Kelvin<=0')
  rho=650+128*C;cp=1450+2736*C/(1+C);S=rho*cp;Sc=128*cp+rho*2736/(1+C)**2
  k=.21+.38*C/(1+C);kc=.38/(1+C)**2;A=.0024*np.exp(-3850/(T+273.15));At=A*3850/(T+273.15)**2
  U=primitive(C);gt=self.G@(T-T[-1]);gu=self.G@(U-U[-1]);qt=k*gt;qm=A*gu
  ht=qt+self.x*(self.G@(qt-qt[-1]));hm=qm+self.x*(self.G@(qm-qm[-1]))
  fT=self.fac*ht/S;fC=self.fac*hm
  if not jac:return np.column_stack([fT[:-1],fC[:-1]]).ravel()
  BU=np.exp(-.45/C)[:,None]*BC
  gtj=self.G@(BT-BT[-1]);guj=self.G@(BU-BU[-1])
  qtj=k[:,None]*gtj+(kc*gt)[:,None]*BC
  qmj=A[:,None]*guj+(At*gu)[:,None]*BT
  htj=qtj+self.x[:,None]*(self.G@(qtj-qtj[-1]));hmj=qmj+self.x[:,None]*(self.G@(qmj-qmj[-1]))
  J=np.empty((2*self.n,2*self.n));J[::2]= (self.fac*htj/S[:,None]-(fT*Sc/S)[:,None]*BC)[:-1];J[1::2]=(self.fac*hmj)[:-1]
  return J
 def profiles(self,u):
  T,C=u.T;U=primitive(C);v=self.I@U;guess=self.I@C
  return np.column_stack([self.I@T,inverse_primitive(v,guess)])
 def extrema(self,u):
  # Interpolate U, not steep C. All real stationary points of the polynomial.
  U=primitive(u[:,1]);co=dct(U[::-1],type=1)/self.n;co[[0,-1]]*=.5
  p=Chebyshev(co);roots=p.deriv().roots();roots=roots[np.isreal(roots)].real;roots=roots[(roots>-1)&(roots<1)]
  pts=np.r_[-1,roots,1];vals=p(pts);idx=np.argmax(vals)
  maxi=inverse_primitive(np.array([vals[idx]]),np.array([u[:,1].max()]))[0]
  return {'max_C':float(maxi),'x_max':float((pts[idx]+1)/2),'stationary_points':len(roots),'minimum_U':float(vals.min())}

def run(n,out,input_file,method='Radau',rtol=2e-11,atol_C=2e-13,max_step=120.):
 if out.with_suffix('.json').exists():raise FileExistsError(out)
 op=Reference(n);data=load_input(input_file);env=Environment(data,'mean');y=np.tile([28.,2.55],n)
 cuts=np.r_[data[:,0],np.arange(21600,360001,21600.)]
 times=[0.];vals=[np.tile([28.,2.55],(21,1))];snapt=[0.];snaps=[np.tile([28.,2.55],(n+1,1))];ext=[]
 nextt=60.;tic=time.perf_counter();nsteps=nfev=nlu=0;event=None;max_dense_overshoot=0.;globalmin=2.55
 for a,b in zip(cuts[:-1],cuts[1:]):
  e=env.segment(a,b)
  fun=lambda t,y:op.evaluate(t,y,e)
  def ev(t,y):return np.max(y[1::2])-.15
  ev.terminal=False;ev.direction=-1
  sol=solve_ivp(fun,(a,b),y,method=method,jac=lambda t,y:op.evaluate(t,y,e,True),rtol=rtol,atol=np.tile([2e-11,atol_C],n),max_step=min(max_step,30.) if a<14400 else max_step,dense_output=True,events=ev)
  if not sol.success:raise RuntimeError(sol.message)
  nsteps+=len(sol.t)-1;nfev+=sol.nfev;nlu+=sol.nlu;end=b
  if len(sol.t_events[0]):
   root=float(sol.t_events[0][0]);end=root+.01;yr=sol.sol(root);ur=op.surface(root,yr.reshape(n,2),e)
   event={'critical_s':root,'critical_h':root/3600,'t_minus_s':root-.01,'max_minus':float(np.max(sol.sol(root-.01)[1::2])),'t_plus_s':end,'max_plus':float(np.max(sol.sol(end)[1::2])),'slope':float(fun(root,yr)[1::2][np.argmax(yr[1::2])]),'polynomial_extrema':op.extrema(ur)}
  ts=np.arange(nextt,end+1e-8,60)
  if event:ts=np.r_[ts,end]
  for t in ts:
   uu=op.surface(t,sol.sol(t).reshape(n,2),e);times.append(float(t));vals.append(op.profiles(uu))
   denseU=op.denseI@primitive(uu[:,1]);max_dense_overshoot=max(max_dense_overshoot,float(denseU.max()-primitive(uu[:,1]).max()));globalmin=min(globalmin,float(uu[:,1].min()))
  if len(ts):nextt=60*(np.floor(ts[-1]/60)+1)
  if b in (1800,3600,5400,7200,9000,10800,14400,21600) or b%21600==0 or event:
   u=op.surface(end,sol.sol(end).reshape(n,2),e);snapt.append(float(end));snaps.append(u);ext.append({'t_s':float(end),**op.extrema(u)})
  y=sol.sol(end)
  if b%21600==0 or event:print(n,b/3600,float(y[1::2].max()),time.perf_counter()-tic,flush=True)
  if event:break
 stats={'n':n,'method':method,'rtol':rtol,'atol_C':atol_C,'max_step':max_step,'event':event,'elapsed_s':time.perf_counter()-tic,'nsteps':nsteps,'nfev':nfev,'nlu':nlu,'min_scalar_boundary_derivative':op.min_boundary_derivative,'max_boundary_residual':op.max_bc_residual,'max_boundary_newton':op.max_newton,'boundary_fallbacks':op.boundary_fallbacks,'dense_U_overshoot':max_dense_overshoot,'global_C_min':globalmin,'extrema':ext}
 out.parent.mkdir(parents=True,exist_ok=True);out.with_suffix('.json').write_text(json.dumps(stats,indent=2));np.savez_compressed(out.with_suffix('.npz'),time_s=np.array(times),sample_TC=np.array(vals),x=op.x,snapshot_time_s=np.array(snapt),snapshots=np.array(snaps))
 print(json.dumps(stats),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=80);p.add_argument('--out',type=Path,required=True);p.add_argument('--method',default='Radau');p.add_argument('--rtol',type=float,default=2e-11);p.add_argument('--input',type=Path,default=Path(__file__).resolve().parents[1]/'inputs/attachment1.xlsx');a=p.parse_args();run(a.n,a.out,a.input,a.method,a.rtol)
