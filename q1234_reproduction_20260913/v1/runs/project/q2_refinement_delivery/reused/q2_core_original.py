"""Q2 effective coupled radial transport. All dimensional quantities are SI.

The temperature storage is S(C)*dT/dt, NOT d(S*T)/dt. The supplied density
is used only inside the empirical effective heat capacity. It is not asserted
to be the actual mass density of a fixed, motionless two-component mixture.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, platform, time, zipfile
from dataclasses import dataclass, asdict
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix

NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

def read_workbook_values(path):
    """Read actual OOXML bytes, including shared strings, without editing input."""
    sheets={}
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None: raise ValueError('Damaged XLSX ZIP')
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS):
                shared.append(''.join(t.text or '' for t in si.iter('{'+NS['s']+'}t')))
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        rel={r.attrib['Id']:r.attrib['Target'] for r in rels}
        for sh in ET.fromstring(z.read('xl/workbook.xml')).find('s:sheets',NS):
            rid=sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            target=rel[rid];target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            vals={}; styles={};types={}
            for c in ET.fromstring(z.read(target)).findall('.//s:sheetData/s:row/s:c',NS):
                addr=c.attrib['r']; typ=c.attrib.get('t','n'); v=c.find('s:v',NS)
                if c.find('s:f',NS) is not None: raise ValueError('Formula in numeric source/result')
                if typ=='s': value=shared[int(v.text)]
                elif typ=='inlineStr': value=''.join(t.text or '' for t in c.findall('.//s:t',NS))
                elif v is None:value=None
                elif typ in ('str','e'):value=v.text
                else:value=float(v.text)
                vals[addr]=value;styles[addr]=int(c.attrib.get('s',0));types[addr]=typ
            sheets[sh.attrib['name']]={'values':vals,'styles':styles,'types':types}
    return sheets

def file_hashes(path):
    b=Path(path).read_bytes()
    return {'sha256':hashlib.sha256(b).hexdigest(),
            'git_blob_sha1':hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest(),
            'bytes':len(b)}

def load_environment(path):
    sheets=read_workbook_values(path)
    if len(sheets)!=1: raise ValueError('Unexpected environmental sheets')
    sh=next(iter(sheets.values()))['values']
    if [sh.get(f'{c}1') for c in 'ABC']!=['时间','温度','水分浓度']:
        raise ValueError('Unrecognized attachment 1 headings')
    last=max(int(a[1:]) for a in sh if a.startswith('A'))
    a=np.array([[sh.get(f'{c}{i}',np.nan) for c in 'ABC'] for i in range(2,last+1)],float)
    if not np.isfinite(a).all() or not np.all(np.diff(a[:,0])>0):
        raise ValueError('Missing, nonnumeric, duplicate, or reversed input')
    if np.any(a[:,2]<0) or a[0,0]!=0:raise ValueError('Invalid initial time/concentration')
    audit={'file':Path(path).name,**file_hashes(path),'sheets':list(sheets),'count':len(a),
           'headings':['时间','温度','水分浓度'],'units':['s','degC','kg/kg (environmental basis unspecified)'],
           'time_range_s':[float(a[0,0]),float(a[-1,0])],
           'time_steps_s':np.unique(np.diff(a[:,0])).tolist(),'missing_values':0,'duplicate_times':0,
           'first_3h_records':int(np.sum(a[:,0]<=10800)),
           'temperature_minmax': [float(a[:,1].min()),float(a[:,1].max())],
           'environmental_moisture_minmax':[float(a[:,2].min()),float(a[:,2].max())],
           'internal_response_observations':False}
    return a,audit

class Environment:
    def __init__(self,data,kind='linear',extension='error'):
        self.data=np.asarray(data,float);self.knots=self.data[:,0]
        if kind not in ('linear','pchip'):raise ValueError(kind)
        if extension not in ('error','hold_last'):raise ValueError(extension)
        self.kind,self.extension=kind,extension
        self.interp=PchipInterpolator(self.knots,self.data[:,1:],axis=0,extrapolate=False) if kind=='pchip' else None
    def __call__(self,t):
        t=np.asarray(t,float)
        if np.any(t<0) or (self.extension=='error' and np.any(t>self.knots[-1]+1e-9)):
            raise ValueError('No implicit environmental extrapolation; explicitly select hold_last for scenarios')
        if self.kind=='pchip':return self.interp(np.minimum(t,self.knots[-1]))
        return np.stack([np.interp(t,self.knots,self.data[:,j]) for j in (1,2)],axis=-1)
    def derivative(self,t):
        t=np.asarray(t,float)
        if self.kind=='pchip':return np.where((t<=self.knots[-1])[...,None],self.interp.derivative()(np.minimum(t,self.knots[-1])),0.)
        ix=np.clip(np.searchsorted(self.knots,t,side='right')-1,0,len(self.knots)-2)
        d=(self.data[ix+1,1:]-self.data[ix,1:])/(self.knots[ix+1]-self.knots[ix])[...,None]
        return np.where((t<self.knots[-1])[...,None],d,0.)

@dataclass(frozen=True)
class Parameters:
    R:float=0.02
    L:float=0.25
    T0:float=28.
    C0:float=2.55
    hT:float=25.
    hm:float=8e-7
P=Parameters()

def properties(T,C,derivatives=False,freeze=False):
    """T in degC, C in kg water/kg dry solid. No silent clipping."""
    T,C=np.broadcast_arrays(np.asarray(T),np.asarray(C))
    if np.any(np.real(C)<=0) or np.any(np.real(T+273.15)<=0):raise FloatingPointError('Nonpositive C or absolute temperature')
    if freeze:T=T*0+28.;C=C*0+2.55
    rho=650.+128.*C;cp=1450.+2736.*C/(1.+C)
    k=.21+.38*C/(1.+C);S=rho*cp
    D=.0024*np.exp(-.45/C-3850./(T+273.15))
    if not derivatives:return S,k,D
    Sc=128.*cp+rho*2736./(1.+C)**2
    kc=.38/(1.+C)**2;Dt=D*3850./(T+273.15)**2;Dc=D*.45/C**2
    if freeze:Sc=Sc*0;kc=kc*0;Dt=Dt*0;Dc=Dc*0
    return S,k,D,Sc,kc,Dt,Dc

def harmonic_with_derivatives(a,b):
    s=a+b
    return 2*a*b/s,2*(b/s)**2,2*(a/s)**2

class CoupledFV:
    """Interleaved [T_0,C_0,T_1,C_1,...], vertex/dual-cell conservative fluxes.
    Temperature state offset supports a genuine degC/K invariance test.
    """
    def __init__(self,n,p=P,temperature_offset=0.,freeze=False,fault=None):
        self.n=int(n);self.m=n+1;self.p=p;self.offset=temperature_offset
        self.freeze=freeze;self.fault=fault
        self.r=np.linspace(0,p.R,n+1);self.dr=p.R/n
        self.faces=np.r_[0.,(self.r[1:]+self.r[:-1])/2,p.R]
        self.w=np.diff(self.faces**2)/2;self.v=self.w.sum();self.s=self.faces[1:-1]/self.dr
    def evaluate(self,t,y,env,with_jac=False):
        u=y.reshape(self.m,2);T=u[:,0]-self.offset;C=u[:,1]
        S,k,D,Sc,kc,Dt,Dc=properties(T,C,True,self.freeze)
        if self.fault=='celsius_exponent':
            D=.0024*np.exp(-.45/C-3850./T);Dt=D*3850./T**2;Dc=D*.45/C**2
        if self.fault=='frozen_D':
            D=np.full_like(C,properties(np.array([28.]),np.array([2.55]))[2][0]);Dt*=0;Dc*=0
        af,akl,akr=harmonic_with_derivatives(k[:-1],k[1:])
        if self.fault=='face_alpha':
            af,_,_=harmonic_with_derivatives((k/S)[:-1],(k/S)[1:])
        df,dl,dr=harmonic_with_derivatives(D[:-1],D[1:])
        dT=T[:-1]-T[1:];dC=C[:-1]-C[1:]
        q=self.s*af*dT;j=self.s*df*dC
        f=np.zeros_like(u);f[:-1,0]-=q;f[1:,0]+=q;f[:-1,1]-=j;f[1:,1]+=j
        ambient=env(t)
        qo=self.p.R*self.p.hT*(T[-1]-ambient[0]);jo=self.p.R*self.p.hm*(C[-1]-ambient[1])
        if self.fault=='face_alpha':qo/=S[-1]
        f[-1]-=[qo,jo];f[:,0]/=self.w*(1 if self.fault=='face_alpha' else S);f[:,1]/=self.w
        if self.fault=='product_storage':f[:,0]-=Sc/S*T*f[:,1]
        if not with_jac:return f.ravel()
        if self.fault is not None:raise ValueError('Fault paths deliberately use numerical Jacobian')
        # 2x2 derivative of each outward edge flux w.r.t. left/right state.
        l=np.empty((self.n,2,2));r=np.empty_like(l)
        l[:,0,0]=self.s*af;r[:,0,0]=-self.s*af
        l[:,0,1]=self.s*akl*kc[:-1]*dT;r[:,0,1]=self.s*akr*kc[1:]*dT
        l[:,1,0]=self.s*dl*Dt[:-1]*dC;r[:,1,0]=self.s*dr*Dt[1:]*dC
        l[:,1,1]=self.s*(df+dl*Dc[:-1]*dC);r[:,1,1]=self.s*(-df+dr*Dc[1:]*dC)
        storage=np.stack([self.w*S,self.w],axis=1)
        ix=np.arange(self.n);rows=[];cols=[];vals=[]
        for ro,sign in ((0,-1.),(1,1.)):
            for co,edge in ((0,l),(1,r)):
                for a in range(2):
                    for b in range(2):
                        rows.append(2*(ix+ro)+a);cols.append(2*(ix+co)+b)
                        vals.append(sign*edge[:,a,b]/storage[ix+ro,a])
        idx=np.arange(self.m)
        rows += [2*idx,np.array([2*self.n,2*self.n+1])]
        cols += [2*idx+1,np.array([2*self.n,2*self.n+1])]
        vals += [-f[:,0]*Sc/S,np.array([-self.p.R*self.p.hT/storage[-1,0],-self.p.R*self.p.hm/storage[-1,1]])]
        jac=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(2*self.m,2*self.m)).tocsc()
        return jac
    def sparsity(self):
        from scipy.sparse import diags,kron
        return kron(diags([np.ones(self.n),np.ones(self.m),np.ones(self.n)],[-1,0,1]),np.ones((2,2)),format='csc')

def dense_derivative(poly,t):
    """Differentiate the accepted BDF dense polynomial, not the ODE RHS."""
    t=np.atleast_1d(t)
    x=(t[None,:]-poly.t_shift[:,None])/poly.denom[:,None]
    prod=np.ones(t.size);der=np.zeros(t.size)
    result=np.zeros((poly.D.shape[1],t.size))
    for j in range(poly.order):
        der=der*x[j]+prod/poly.denom[j];prod*=x[j]
        result+=poly.D[j+1,:,None]*der
    return result

def solve_fv(env,n=320,p=P,end=10800.,rtol=2e-11,atol_T=2e-12,atol_C=2e-13,max_step=20.,
             temperature_offset=0.,freeze=False,initial=None,quadrature=False,fault=None,output_step=1.,verbose=False):
    if n%20:raise ValueError('FV n must be divisible by 20 to sample requested radii without interpolation')
    if end>env.knots[-1] and env.extension=='error':raise ValueError('End outside data range')
    op=CoupledFV(n,p,temperature_offset,freeze,fault)
    y=np.empty((n+1,2));y[:,0]=p.T0+temperature_offset;y[:,1]=p.C0
    if initial is not None:y=np.array(initial(op.r),float);y[:,0]+=temperature_offset
    times=np.unique(np.r_[np.arange(0,end+1e-9,output_step),end]);oi=np.arange(21)*(n//20)
    out=np.empty((len(times),21,2));out[0]=y[oi];out[0,:,0]-=temperature_offset
    average=np.empty((len(times),2));average[0]=op.w@y/op.v;average[0,0]-=temperature_offset
    snaps=[y.copy()];st=[0.];stats={'n':n,'nodes':n+1,'dr_m':op.dr,'end_s':end,'rtol':rtol,
      'atol_T':atol_T,'atol_C':atol_C,'max_step_s':max_step,'method':'BDF + exact coupled block sparse Jacobian',
      'environment':env.kind,'extension':env.extension,'temperature_offset':temperature_offset,
      'freeze':freeze,'fault':fault,'nfev':0,'njev':0,'nlu':0,'accepted_steps':0,'parameters':asdict(p)}
    cuts=np.unique(np.r_[0.,env.knots[(env.knots>0)&(env.knots<end)],end])
    tic=time.perf_counter();qint=jint=storeint=0.;balance=[];field_min=np.array([np.inf,np.inf]);field_max=-field_min
    gx,gw=np.polynomial.legendre.leggauss(6)
    atol=np.tile([atol_T,atol_C],n+1)
    for a,b in zip(cuts[:-1],cuts[1:]):
        fun=lambda t,z:op.evaluate(t,z,env)
        jac=(lambda t,z:op.evaluate(t,z,env,True)) if fault is None else None
        sol=solve_ivp(fun,(a,b),y.ravel(),method='BDF',jac=jac,jac_sparsity=op.sparsity() if fault else None,
          rtol=rtol,atol=atol,max_step=max_step,dense_output=True,first_step=min(1e-3,.05*op.dr**2/2e-8,b-a))
        if not sol.success:raise RuntimeError(f'FV failed n={n}, [{a},{b}]: {sol.message}')
        if not np.isfinite(sol.y).all() or np.any(sol.y[1::2]<=0):raise FloatingPointError('Invalid accepted state')
        for key in ('nfev','njev','nlu'):stats[key]+=int(getattr(sol,key))
        stats['accepted_steps']+=len(sol.t)-1
        states=sol.y.reshape(n+1,2,-1)
        field_min=np.minimum(field_min,states.min(axis=(0,2)));field_max=np.maximum(field_max,states.max(axis=(0,2)))
        mask=(times>a)&(times<=b+1e-8)
        v=sol.sol(times[mask]).reshape(n+1,2,-1)
        out[mask]=v[oi].transpose(2,0,1);out[mask,:,0]-=temperature_offset
        average[mask]=np.einsum('i,ijt->tj',op.w,v)/op.v;average[mask,0]-=temperature_offset
        for t in (1800.,3600.,5400.,7200.,9000.,10800.):
            if a<t<=b:snaps.append(sol.sol(t).reshape(n+1,2));st.append(t)
        if quadrature:
            for lo,hi,poly in zip(sol.t[:-1],sol.t[1:],sol.sol.interpolants):
                tq=(lo+hi)/2+(hi-lo)/2*gx;ys=poly(tq).reshape(n+1,2,-1)
                yd=dense_derivative(poly,tq).reshape(n+1,2,-1)
                S=properties(ys[:,0]-temperature_offset,ys[:,1],freeze=freeze)[0]
                boundary=env(tq)
                qdot=p.R*p.hT*(ys[-1,0]-temperature_offset-boundary[:,0])/op.v
                jdot=p.R*p.hm*(ys[-1,1]-boundary[:,1])/op.v
                # Independently differentiate accepted interpolation; no flux cancellation RHS reuse.
                storedot=np.einsum('i,it,it->t',op.w,S,yd[:,0])/op.v
                weights=gw*(hi-lo)/2
                qint+=float(qdot@weights);jint+=float(jdot@weights);storeint+=float(storedot@weights)
            currentC=float(op.w@sol.y[1::2,-1]/op.v)
            balance.append([b,currentC-average[0,1]+jint,storeint+qint,qint,storeint])
        y=sol.y[:,-1].reshape(n+1,2)
        if verbose and int(b)%1800==0:print(f'FV n={n}, t={b:g} s',flush=True)
    stats['elapsed_s']=time.perf_counter()-tic
    field_min[0]-=temperature_offset;field_max[0]-=temperature_offset
    stats['accepted_internal_min']=field_min.tolist();stats['accepted_internal_max']=field_max.tolist()
    stats['environment_knots_restarted']=len(cuts)-1
    stats['quadrature_performed']=quadrature
    if balance:
        ba=np.array(balance);stats['water_balance_max_kgkg']=float(np.max(abs(ba[:,1])))
        stats['effective_heat_balance_max_J_m3']=float(np.max(abs(ba[:,2])))
        stats['net_heat_input_J_m3']=-qint;stats['effective_heat_balance_relative']=float(np.max(abs(ba[:,2]))/max(abs(qint),1.))
    return {'time_s':times,'radius_cm':np.linspace(0,2,21),'temperature_degC':out[:,:,0],
       'moisture_dry_basis':out[:,:,1],'average_temperature_degC':average[:,0],
       'average_moisture_dry_basis':average[:,1],'environment':env(times),'internal_radius_m':op.r,
       'internal_weights_rdr':op.w,'snapshot_times_s':np.array(st),
       'snapshots_temperature_degC':np.array(snaps)[:,:,0]-temperature_offset,
       'snapshots_moisture_dry_basis':np.array(snaps)[:,:,1],'balance_records':np.array(balance),'stats':stats}

def save_solution(path,sol):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**{k:v for k,v in sol.items() if k!='stats'})
    path.with_suffix('.json').write_text(json.dumps(sol['stats'],ensure_ascii=False,indent=2),encoding='utf-8')

def main():
    base=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=base/'inputs/附件1.xlsx');ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n',type=int,default=1280);ap.add_argument('--end',type=float,default=10800.)
    ap.add_argument('--rtol',type=float,default=2e-11);ap.add_argument('--max-step',type=float,default=20.)
    ap.add_argument('--quadrature',action='store_true');ap.add_argument('--extension',choices=['error','hold_last'],default='error')
    args=ap.parse_args();a,audit=load_environment(args.input);args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'input_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    sol=solve_fv(Environment(a,extension=args.extension),n=args.n,end=args.end,rtol=args.rtol,max_step=args.max_step,quadrature=args.quadrature,verbose=True)
    save_solution(args.out/f'fv_n{args.n}.npz',sol);print(json.dumps(sol['stats'],indent=2))
if __name__=='__main__':main()
