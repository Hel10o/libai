"""Q4 material-coordinate Chebyshev collocation, adapted from audited Q3 reference.
Changes: Appendix 4 local coefficients; radius(t); material-coordinate derivation.
Source commit anonymous-source.
Independent global Chebyshev collocation with primitive-gradient moisture flux.
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
from q4_common import load_input,Environment

def primitive(c):
 if np.any(np.real(c)<=0):raise FloatingPointError('C<=0 in reference')
 return c*np.exp(-.30/c)-.30*exp1(.30/c)

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
  delta=(primitive(c)-v)/np.exp(-.30/c)
  cnew=c-delta
  # Globalize Newton without altering the equation or accepted value.
  bad=np.real(cnew)<=0
  while np.any(bad):delta[bad]*=.5;cnew=c-delta;bad=np.real(cnew)<=0
  c=cnew
  if np.max(abs(delta))<3e-14:break
 else:raise RuntimeError('Inverse primitive failed')
 return c

class Reference:
 def __init__(self,n,radius):
  self.radius=radius
  self.n=n;self.x,self.G,self.w=grid(n);self.row=self.G[-1,:-1];self.diag=self.G[-1,-1]
  self.fac=4/.02**2;self.bd=2/.02
  self.BT0=np.zeros((n+1,2*n));self.BC0=np.zeros_like(self.BT0)
  self.BT0[np.arange(n),2*np.arange(n)]=1;self.BC0[np.arange(n),2*np.arange(n)+1]=1
  self.I=interp(self.x,self.w,np.linspace(0,1,21)**2)
  self.denseI=interp(self.x,self.w,np.linspace(0,1,8*n+1))
  self.min_boundary_derivative=float('inf');self.max_bc_residual=0.;self.max_newton=0;self.boundary_fallbacks=0
 def surface(self,t,u,env,jac=False):
  R=self.radius(t);self.fac=4/R**2;self.bd=2/R
  Ti,Ci=u.T;et,ec=env(t);Ui=primitive(Ci);anchor=Ci[-1];ua=Ui[-1]
  gt=self.row@(Ti-et);gu=self.row@(Ui-ua)
  def calc(c):
   k=.12+.20*c/(1+c);kc=.20/(1+c)**2
   den=self.bd*k*self.diag+25
   Ts=et-self.bd*k*gt/den
   Tsc=-self.bd*kc*gt*25/den**2
   A=.00042*np.exp(-3850/(Ts+273.15));At=A*3850/(Ts+273.15)**2
   ux=gu+self.diag*(primitive(c)-ua)
   f=np.exp(-.30/c)
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
  B=np.array([[self.bd*k*self.diag+25,self.bd*kc*tx],[self.bd*At*ux,self.bd*A*self.diag*np.exp(-.30/c)+8e-7]])
  Q=np.zeros((2,2*self.n));Q[0,::2]=self.bd*k*self.row;Q[1,1::2]=self.bd*A*self.row*np.exp(-.30/Ci)
  sens=-np.linalg.solve(B,Q);BT=self.BT0.copy();BC=self.BC0.copy();BT[-1]=sens[0];BC[-1]=sens[1]
  return full,BT,BC
 def evaluate(self,t,y,env,jac=False):
  if jac:full,BT,BC=self.surface(t,y.reshape(self.n,2),env,True)
  else:full=self.surface(t,y.reshape(self.n,2),env)
  T,C=full.T
  if np.any(T<=-273.15):raise FloatingPointError('Kelvin<=0')
  rho=760+90*C;cp=1850+2150*C/(1+C);S=rho*cp;Sc=90*cp+rho*2150/(1+C)**2
  k=.12+.20*C/(1+C);kc=.20/(1+C)**2;A=.00042*np.exp(-3850/(T+273.15));At=A*3850/(T+273.15)**2
  U=primitive(C);gt=self.G@(T-T[-1]);gu=self.G@(U-U[-1]);qt=k*gt;qm=A*gu
  ht=qt+self.x*(self.G@(qt-qt[-1]));hm=qm+self.x*(self.G@(qm-qm[-1]))
  fT=self.fac*ht/S;fC=self.fac*hm
  if not jac:return np.column_stack([fT[:-1],fC[:-1]]).ravel()
  BU=np.exp(-.30/C)[:,None]*BC
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
