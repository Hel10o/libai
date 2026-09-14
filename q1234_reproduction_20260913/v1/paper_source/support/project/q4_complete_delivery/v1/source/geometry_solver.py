"""Independent vertex/dual-volume axisymmetric validation of Q1/Q23/Q4.

All states remain coupled from t=0. No isothermal projection is used. The
reference grid follows homogeneous radial material contraction, while length
is fixed. Reference weights are constant; radial diffusion scales as q^-2,
side exchange as q^-1, axial diffusion/end exchange remain unchanged.
"""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix
from scipy.special import exp1
from openpyxl import load_workbook

ROOT=Path(__file__).resolve().parents[3]
DEST=Path(__file__).resolve().parents[1]/'validation'/'geometry'

def data_xlsx(path, width):
    wb=load_workbook(path,read_only=True,data_only=True)
    rows=list(wb.active.values); wb.close()
    a=np.asarray([row[:width] for row in rows[1:] if row[0] is not None],float)
    if not np.isfinite(a).all() or not np.all(np.diff(a[:,0])>0):
        raise ValueError('Invalid input records')
    return a

class Inputs:
    def __init__(self,input_dir=None):
        packaged=Path(input_dir).resolve() if input_dir else Path(__file__).resolve().parents[1]/'inputs'
        if input_dir is not None or (packaged/'attachment1.xlsx').exists():
            env_path=packaged/'attachment1.xlsx';radius_path=packaged/'attachment2.xlsx'
        else:
            env_path=ROOT/'A题/附件/附件1.xlsx';radius_path=ROOT/'A题/附件/附件2.xlsx'
        self.env=data_xlsx(env_path,3)
        self.rad=data_xlsx(radius_path,2)
        self.rad[:,1]*=.01
        self.future=self.env[self.env[:,0]>=10800,1:].mean(axis=0)
    def ambient(self,t,future=False):
        if future:return self.future
        return np.array([np.interp(t,self.env[:,0],self.env[:,j]) for j in (1,2)])
    def radius(self,t):
        if t>self.rad[-1,0]+1e-6:raise ValueError('Radius input ends at 72 h')
        return float(np.interp(t,self.rad[:,0],self.rad[:,1]))

def properties(T,C,case):
    if np.any(C<=0) or np.any(T<=-273.15):raise FloatingPointError('Nonpositive material state')
    if case=='q1':
        S=np.full_like(C,820.*2600.); k=np.full_like(C,.36)
        D=7e-9*np.exp(-.89/C)
        return S,k,D,np.zeros_like(C),np.zeros_like(C),np.zeros_like(C),D*.89/C**2
    if case=='q23':a,b,c,d,e,f,A,B=650.,128.,1450.,2736.,.21,.38,.0024,.45
    elif case=='q4':a,b,c,d,e,f,A,B=760.,90.,1850.,2150.,.12,.20,.00042,.30
    else:raise ValueError(case)
    rho=a+b*C;cp=c+d*C/(1+C);S=rho*cp;k=e+f*C/(1+C)
    D=A*np.exp(-B/C-3850/(T+273.15))
    return S,k,D,b*cp+rho*d/(1+C)**2,f/(1+C)**2,D*3850/(T+273.15)**2,D*B/C**2

def primitive_difference(a,b,B):
    mid=(a+b)/2;delta=a-b;close=np.abs(delta)<.002*mid;value=np.empty_like(a)
    if np.any(~close):
        x=a[~close];y=b[~close]
        value[~close]=x*np.exp(-B/x)-B*exp1(B/x)-y*np.exp(-B/y)+B*exp1(B/y)
    if np.any(close):
        m=mid[close];d=delta[close];v=np.zeros_like(m)
        for x,w in ((-.8611363115940526,.3478548451374539),(-.3399810435848563,.6521451548625461),(.3399810435848563,.6521451548625461),(.8611363115940526,.3478548451374539)):
            v+=w*np.exp(-B/(m+x*d/2))
        value[close]=v*d/2
    return value

class GeometryFV:
    def __init__(self,nr,nz,case,inputs,grading=2.,zgrading=2.,end_transfer=True,shrink=True,flux='integral'):
        self.nr,self.nz,self.case,self.inputs=nr,nz,case,inputs
        self.shrink=shrink and case=='q4'
        self.flux=flux
        self.r=.02*(1-(1-np.linspace(0,1,nr+1))**grading)
        rf=np.r_[0,(self.r[:-1]+self.r[1:])/2,.02];wr=np.diff(rf**2)/2
        if nz:
            self.z=.125*(1-(1-np.linspace(0,1,nz+1))**zgrading)
            zf=np.r_[0,(self.z[:-1]+self.z[1:])/2,.125];wz=np.diff(zf)
        else:self.z=np.array([0.]);wz=np.ones(1)
        self.shape=(nz+1,nr+1);self.m=(nr+1)*(nz+1)
        self.w=(wz[:,None]*wr).ravel();idx=np.arange(self.m).reshape(self.shape)
        i=[idx[:,:-1].ravel()];j=[idx[:,1:].ravel()]
        g=[(wz[:,None]*rf[None,1:-1]/np.diff(self.r)).ravel()]
        self.radial_edges=len(i[0])
        if nz:
            i.append(idx[:-1].ravel());j.append(idx[1:].ravel())
            g.append((wr[None,:]/np.diff(self.z)[:,None]).ravel())
        self.i=np.concatenate(i);self.j=np.concatenate(j);self.g=np.concatenate(g)
        side=np.zeros(self.shape);end=np.zeros(self.shape)
        side[:,-1]=.02*wz
        if nz and end_transfer:end[-1,:]=wr
        self.side=side.ravel();self.end=end.ravel()
        rows=[];cols=[]
        for dest in (self.i,self.j):
            for src in (self.i,self.j):
                for a in range(2):
                    for b in range(2):rows.append(2*dest+a);cols.append(2*src+b)
        nodes=np.arange(self.m)
        self.jrows=np.r_[np.concatenate(rows),2*nodes,2*nodes,2*nodes+1]
        self.jcols=np.r_[np.concatenate(cols),2*nodes+1,2*nodes,2*nodes+1]
    def rhs(self,t,y,ambient,jac=False):
        q=self.inputs.radius(t)/.02 if self.shrink else 1.
        g=self.g.copy();g[:self.radial_edges]/=q*q
        area=self.side/q+self.end
        U=y.reshape(self.m,2);T,C=U.T;i,j=self.i,self.j
        S,k,D,Sc,kc,Dt,Dc=properties(T,C,self.case)
        ks=k[i]+k[j];kf=2*k[i]*k[j]/ks;kl=2*(k[j]/ks)**2;kr=2*(k[i]/ks)**2
        ds=D[i]+D[j];df=2*D[i]*D[j]/ds;dl=2*(D[j]/ds)**2;dr=2*(D[i]/ds)**2
        dT=T[i]-T[j];dC=C[i]-C[j];fT=g*kf*dT
        if self.flux=='integral':
            B=.89 if self.case=='q1' else .45 if self.case=='q23' else .30
            A=7e-9 if self.case=='q1' else (.0024 if self.case=='q23' else .00042)*np.exp(-3850/((T[i]+T[j])/2+273.15))
            fC=g*A*primitive_difference(C[i],C[j],B)
            mci=g*A*np.exp(-B/C[i]);mcj=-g*A*np.exp(-B/C[j])
            mti=np.zeros_like(fC) if self.case=='q1' else fC*1925/((T[i]+T[j])/2+273.15)**2
            mtj=mti
        else:
            fC=g*df*dC
            mci=g*(df+dl*Dc[i]*dC);mcj=g*(-df+dr*Dc[j]*dC)
            mti=g*dl*Dt[i]*dC;mtj=g*dr*Dt[j]*dC
        f=np.column_stack([np.bincount(j,weights=fT,minlength=self.m)-np.bincount(i,weights=fT,minlength=self.m),
                           np.bincount(j,weights=fC,minlength=self.m)-np.bincount(i,weights=fC,minlength=self.m)])
        et,ec=ambient(t);f[:,0]-=area*25.*(T-et);f[:,1]-=area*8e-7*(C-ec)
        storage=self.w[:,None]*np.column_stack([S,np.ones(self.m)]);f/=storage
        if not jac:return f.ravel()
        left=np.empty((len(i),2,2));right=np.empty_like(left)
        left[:,0,0]=g*kf;right[:,0,0]=-g*kf
        left[:,0,1]=g*kl*kc[i]*dT;right[:,0,1]=g*kr*kc[j]*dT
        left[:,1,0]=mti;right[:,1,0]=mtj
        left[:,1,1]=mci;right[:,1,1]=mcj
        vals=[]
        for dest,sgn in ((i,-1.),(j,1.)):
            for v in (left,right):
                for a in range(2):
                    for b in range(2):vals.append(sgn*v[:,a,b]/storage[dest,a])
        vals.extend([-f[:,0]*Sc/S,-area*25./storage[:,0],-area*8e-7/self.w])
        return coo_matrix((np.concatenate(vals),(self.jrows,self.jcols)),shape=(2*self.m,2*self.m)).tocsc()
    def observe(self,t,y):
        u=y.reshape(self.shape+(2,));C=u[:,:,1];ix=np.unravel_index(np.argmax(C),C.shape)
        radius=self.inputs.radius(t) if self.shrink else .02
        # Material-coordinate samples, exactly same locations for paired 1D/2D.
        target=np.linspace(0,.02,21)
        mid=PchipInterpolator(self.r**2,u[0],axis=0)(target**2)
        end=PchipInterpolator(self.r**2,u[-1],axis=0)(target**2)
        return mid,end,float(C[ix]),np.array([self.r[ix[1]]*radius/.02,self.z[ix[0]]]),radius

def run(case,nr,nz,end=None,rtol=2e-9,max_step=240.,shrink=True,end_transfer=True,flux='integral',out_dir=None,input_dir=None):
    source_hash_at_start=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest=Path(out_dir).resolve() if out_dir else DEST
    inp=Inputs(input_dir);op=GeometryFV(nr,nz,case,inp,shrink=shrink,end_transfer=end_transfer,flux=flux)
    end=float(end if end is not None else (1800 if case=='q1' else 220000 if case=='q23' else 259200))
    label=f'{case}_{nr}x{nz}'+('_fixed' if case=='q4' and not shrink else '')+('_sealed_ends' if nz and not end_transfer else '')+('_integral' if flux=='integral' else '')
    if rtol!=2e-9:label+=f'_rtol{rtol:g}'
    dest.mkdir(parents=True,exist_ok=True)
    if any((dest/(label+suffix)).exists() for suffix in ('.json','.npz')):raise FileExistsError(label)
    start=time.perf_counter();y=np.tile([28.,2.55],op.m)
    sample_times=np.unique(np.r_[np.arange(0,end+1,60.),end,[100,300,900,1200,1500]])
    sample_times=sample_times[(sample_times>=0)&(sample_times<=end)]
    snapshots_at=np.unique(np.r_[0,100,300,600,900,1200,1500,1800,3600,5400,7200,9000,10800,14400,np.arange(21600,end,21600)])
    cuts=np.unique(np.r_[0,inp.env[:,0][inp.env[:,0]<end],inp.rad[:,0][inp.rad[:,0]<end] if op.shrink else [],end])
    recorded=[];mid=[];edge=[];cmax=[];locations=[];radii=[];snapshots=[y.reshape(op.shape+(2,)).copy()];st=[0.]
    stats={'case':case,'nr':nr,'nz':nz,'grading':2.,'zgrading':2.,'shrink':op.shrink,'end_transfer':end_transfer,
           'rtol':rtol,'atol_T':1e-9,'atol_C':1e-11,'max_step':max_step,'nfev':0,'njev':0,'nlu':0,
           'flux':flux+' local coefficients, vertex dual finite volumes','temperature':'coupled throughout; no isothermal projection',
           'environment_future':inp.future.tolist(),'source_sha':'anonymous-source'}
    def record(t,Y):
        m,e,c,l,R=op.observe(t,Y);recorded.append(t);mid.append(m);edge.append(e);cmax.append(c);locations.append(l);radii.append(R)
    record(0.,y)
    def event(t,Y):return np.max(Y[1::2])-.15
    event.direction=-1;event.terminal=True
    root=None;event_states=[];event_times=[]
    for a,b in zip(cuts[:-1],cuts[1:]):
        ambient=(lambda t: inp.ambient(t,True)) if a>=14400 else inp.ambient
        sol=solve_ivp(lambda t,Y:op.rhs(t,Y,ambient),(a,b),y,method='BDF',
            jac=lambda t,Y:op.rhs(t,Y,ambient,True),rtol=rtol,atol=np.tile([1e-9,1e-11],op.m),
            max_step=min(max_step,30.) if a<14400 else max_step,dense_output=True,
            events=event if case!='q1' else None)
        if not sol.success:raise RuntimeError(sol.message)
        for key in ('nfev','njev','nlu'):stats[key]+=int(getattr(sol,key))
        stop=float(sol.t[-1]);ts=sample_times[(sample_times>a)&(sample_times<=stop)]
        if len(ts):
            vals=sol.sol(ts)
            for t,Y in zip(ts,vals.T):record(float(t),Y)
        sts=snapshots_at[(snapshots_at>a)&(snapshots_at<=stop)]
        if len(sts):
            for t,Y in zip(sts,sol.sol(sts).T):st.append(float(t));snapshots.append(Y.reshape(op.shape+(2,)).copy())
        y=sol.y[:,-1]
        if sol.t_events is not None and len(sol.t_events[0]):
            root=float(sol.t_events[0][0]);record(root,y);st.append(root);snapshots.append(y.reshape(op.shape+(2,)).copy())
            before=max(float(a),root-1.);after=root+1.
            ext=solve_ivp(lambda t,Y:op.rhs(t,Y,ambient),(root,after),y,method='BDF',
                jac=lambda t,Y:op.rhs(t,Y,ambient,True),rtol=rtol,atol=np.tile([1e-9,1e-11],op.m),dense_output=True)
            if not ext.success:raise RuntimeError(ext.message)
            event_times=[before,root,after];event_states=[sol.sol(before).reshape(op.shape+(2,)),y.reshape(op.shape+(2,)),ext.y[:,-1].reshape(op.shape+(2,))]
            break
        if b%21600==0 or b in (1800,10800,14400):print(label,'t=',b,'max C=',float(np.max(y[1::2])),'elapsed=',round(time.perf_counter()-start,2),flush=True)
    stats['elapsed_s']=time.perf_counter()-start;stats['event_root_s']=root
    stats['end_s']=float(recorded[-1]);stats['end_max_C']=float(cmax[-1]);stats['end_max_location_r_z_m']=locations[-1].tolist()
    if root is not None:stats['event_bracket_Cmax']=[float(np.max(v[:,:,1])) for v in event_states]
    stats['source_sha256']=source_hash_at_start
    np.savez_compressed(dest/(label+'.npz'),time_s=recorded,mid=mid,end=edge,Cmax=cmax,Cmax_location_r_z_m=locations,
        radius_m=radii,reference_r_m=op.r,z_m=op.z,weights=op.w,snapshot_times_s=st,snapshots=snapshots,
        event_times_s=event_times,event_states=event_states)
    (dest/(label+'.json')).write_text(json.dumps(stats,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(stats,ensure_ascii=False),flush=True)
    return stats

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',choices=['q1','q23','q4'],required=True)
    p.add_argument('--nr',type=int,required=True);p.add_argument('--nz',type=int,default=0)
    p.add_argument('--end',type=float);p.add_argument('--rtol',type=float,default=2e-9)
    p.add_argument('--max-step',type=float,default=240.);p.add_argument('--fixed',action='store_true');p.add_argument('--sealed-ends',action='store_true')
    p.add_argument('--flux',choices=['integral','harmonic'],default='integral')
    p.add_argument('--out-dir',type=Path,help='New output location; existing run files are never overwritten.')
    p.add_argument('--input-dir',type=Path,help='Directory containing attachment1.xlsx and attachment2.xlsx.')
    a=p.parse_args();run(a.case,a.nr,a.nz,a.end,a.rtol,a.max_step,not a.fixed,not a.sealed_ends,a.flux,a.out_dir,a.input_dir)
