"""Q3: independently assembled conservative graph FV on graded cylindrical meshes.
SI units; T in Celsius except Arrhenius; C is dry-basis. No clipping or D floor.
Both actual early boundary knots and the future scenario are explicit.
"""
from __future__ import annotations
import argparse, json, time, sys, hashlib, platform, zipfile
from pathlib import Path
from dataclasses import dataclass, asdict
import xml.etree.ElementTree as ET
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import exp1
from scipy.optimize import brentq
from scipy.sparse import coo_matrix

NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

def read_xlsx(path):
    """Read physical OOXML (not preview text), retaining numeric/string identity."""
    with zipfile.ZipFile(path) as z:
        if z.testzip(): raise ValueError('Corrupt xlsx')
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for v in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS):
                strings.append(''.join(t.text or '' for t in v.findall('.//s:t',NS)))
        rel={v.attrib['Id']:v.attrib['Target'] for v in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        out={}
        for sh in ET.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',NS):
            rid=sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            target=rel[rid]; target=target[1:] if target.startswith('/') else 'xl/'+target
            cells={}
            for c in ET.fromstring(z.read(target)).findall('.//s:sheetData/s:row/s:c',NS):
                kind=c.attrib.get('t','n');v=c.find('s:v',NS)
                if c.find('s:f',NS) is not None: raise ValueError('Unexpected formula')
                if kind=='s': value=strings[int(v.text)]
                elif kind=='inlineStr': value=''.join(t.text or '' for t in c.findall('.//s:t',NS))
                elif v is None: value=None
                elif kind in ('str','e'): value=v.text
                else: value=float(v.text)
                cells[c.attrib['r']]={'value':value,'type':kind,'style':int(c.attrib.get('s','0'))}
            out[sh.attrib['name']]=cells
        return out

def load_input(path):
    cells=read_xlsx(path)['Sheet1']
    assert [cells[f'{c}1']['value'] for c in 'ABC']==['时间','温度','水分浓度']
    a=np.array([[cells[f'{c}{i}']['value'] for c in 'ABC'] for i in range(2,243)],float)
    assert a.shape==(241,3) and np.isfinite(a).all()
    assert np.array_equal(a[:,0],np.arange(241)*60.)
    return a

class Environment:
    def __init__(self,data,scenario='mean',kind='linear'):
        self.data=np.asarray(data,float);self.scenario=scenario;self.kind=kind
        if scenario=='mean':self.future=self.data[self.data[:,0]>=10800,1:].mean(axis=0)
        elif scenario=='last':self.future=self.data[-1,1:].copy()
        elif scenario=='nominal':self.future=np.array([50.,.05])
        else:raise ValueError(scenario)
        self.interp=PchipInterpolator(self.data[:,0],self.data[:,1:],axis=0) if kind=='pchip' else None
    def __call__(self,t):
        t=np.asarray(t);v=(self.interp(t) if self.interp is not None else np.stack([np.interp(t,self.data[:,0],self.data[:,k]) for k in (1,2)],axis=-1))
        return np.where((t>14400)[...,None],self.future,v)
    def segment(self,a,b):
        # One-sided boundary at t=14400: no spurious jump in the pre-4h segment.
        if a>=14400:return lambda t: np.broadcast_to(self.future,np.shape(t)+(2,))
        return self

def props(T,C):
    if np.any(C<=0) or np.any(T<=-273.15):raise FloatingPointError('Nonpositive accepted/trial C or Kelvin temperature')
    rho=650.+128*C;cp=1450.+2736*C/(1+C)
    S=rho*cp;k=.21+.38*C/(1+C)
    Sc=128*cp+rho*2736/(1+C)**2;kc=.38/(1+C)**2
    A=.0024*np.exp(-3850/(T+273.15));D=A*np.exp(-.45/C)
    return S,k,D,Sc,kc

def F(C):
    C=np.asarray(C)
    if np.any(C<=0):raise FloatingPointError('F requires C>0')
    return C*np.exp(-.45/C)-.45*exp1(.45/C)

def Fdiff(a,b):
    """Stable primitive difference. 4-point Gaussian integration when close.
    The switch is numerical evaluation only; no modification of diffusivity.
    """
    a,b=np.broadcast_arrays(a,b);mid=(a+b)/2;delta=a-b
    close=np.abs(delta)<.002*mid
    ans=np.empty_like(a)
    if np.any(~close):ans[~close]=F(a[~close])-F(b[~close])
    if np.any(close):
        m=mid[close];d=delta[close]
        v=np.zeros_like(m)
        for x,w in ((-.8611363115940526,.3478548451374539),(-.3399810435848563,.6521451548625461),(.3399810435848563,.6521451548625461),(.8611363115940526,.3478548451374539)):
            v+=w*np.exp(-.45/(m+d*x/2))
        ans[close]=d*v/2
    return ans

@dataclass(frozen=True)
class Config:
    n:int=256
    nz:int=0
    grading:float=2.
    zgrading:float=2.
    flux:str='integral'
    scenario:str='mean'
    interpolation:str='linear'
    hT:float=25.
    hm:float=8e-7
    end_factor:float=1.
    rtol:float=2e-9
    atol_T:float=2e-10
    atol_C:float=2e-12
    max_step:float=240.
    early_step:float=30.
    method:str='BDF'
    end_s:float=360000.
    isothermal_after:float=0.
    quadrature:bool=False
    full_samples:bool=False
    execution_time_s:float=0.

class Operator:
    def __init__(self,cfg:Config):
        self.cfg=cfg;self.R=.02;self.Lhalf=.125
        self.r=self.R*(1-(1-np.linspace(0,1,cfg.n+1))**cfg.grading)
        rf=np.r_[0,(self.r[:-1]+self.r[1:])/2,self.R];wr=np.diff(rf**2)/2
        if cfg.nz:
            self.z=self.Lhalf*(1-(1-np.linspace(0,1,cfg.nz+1))**cfg.zgrading)
            zf=np.r_[0,(self.z[:-1]+self.z[1:])/2,self.Lhalf];wz=np.diff(zf)
        else:self.z=np.array([0.]);wz=np.ones(1)
        self.shape=(len(self.r),len(self.z)); self.size=np.prod(self.shape).item()
        ids=np.arange(self.size).reshape(self.shape)
        self.w=(wr[:,None]*wz).ravel();self.vol=self.w.sum()
        l=[ids[:-1].ravel()];r=[ids[1:].ravel()]
        ge=[((rf[1:-1]/np.diff(self.r))[:,None]*wz).ravel()]
        if cfg.nz:
            l.append(ids[:,:-1].ravel());r.append(ids[:,1:].ravel())
            ge.append((wr[:,None]/np.diff(self.z)).ravel())
        self.left=np.concatenate(l);self.right=np.concatenate(r);self.geo=np.concatenate(ge)
        self.bound=np.zeros(self.shape);self.bound[-1]+=self.R*wz
        if cfg.nz:self.bound[:,-1]+=cfg.end_factor*wr
        self.bound=self.bound.ravel()
        # Fixed sparse block pattern, row divided by its own storage.
        rows=[];cols=[]; edgeid=[]; signs=[]
        for dest,sgn in ((self.left,-1.),(self.right,1.)):
            for src in (self.left,self.right):
                for a in range(2):
                    for b in range(2):
                        rows.append(2*dest+a);cols.append(2*src+b)
        self.jrows=np.r_[np.concatenate(rows),2*np.arange(self.size),2*np.arange(self.size),2*np.arange(self.size)+1]
        self.jcols=np.r_[np.concatenate(cols),2*np.arange(self.size)+1,2*np.arange(self.size),2*np.arange(self.size)+1]
    def rhs(self,t,y,env,jac=False):
        q=y.reshape(self.size,2);T=q[:,0];C=q[:,1]
        S,k,D,Sc,kc=props(T,C);l=self.left;r=self.right;g=self.geo
        kl=k[l];kr=k[r];ks=kl+kr;kf=2*kl*kr/ks
        dT=T[l]-T[r];dC=C[l]-C[r]
        tq=g*kf*dT
        if self.cfg.flux=='integral':
            Tf=(T[l]+T[r])/2;Af=.0024*np.exp(-3850/(Tf+273.15));at=3850/(Tf+273.15)**2/2
            mf=g*Af*Fdiff(C[l],C[r])
            mcl=g*Af*np.exp(-.45/C[l]);mcr=-g*Af*np.exp(-.45/C[r])
            mtl=mf*at;mtr=mtl
        else:
            Dl=D[l];Dr=D[r]
            if self.cfg.flux=='harmonic':
                ds=Dl+Dr;df=2*Dl*Dr/ds;fl=2*(Dr/ds)**2;fr=2*(Dl/ds)**2
            elif self.cfg.flux=='arithmetic':df=(Dl+Dr)/2;fl=fr=.5
            else:raise ValueError(self.cfg.flux)
            mf=g*df*dC
            mcl=g*(df+fl*Dl*.45/C[l]**2*dC);mcr=g*(-df+fr*Dr*.45/C[r]**2*dC)
            mtl=g*fl*Dl*3850/(T[l]+273.15)**2*dC;mtr=g*fr*Dr*3850/(T[r]+273.15)**2*dC
        rhs=np.empty_like(q)
        rhs[:,0]=np.bincount(r,weights=tq,minlength=self.size)-np.bincount(l,weights=tq,minlength=self.size)
        rhs[:,1]=np.bincount(r,weights=mf,minlength=self.size)-np.bincount(l,weights=mf,minlength=self.size)
        et,ec=env(t);rhs[:,0]-=self.bound*self.cfg.hT*(T-et);rhs[:,1]-=self.bound*self.cfg.hm*(C-ec)
        storage=self.w[:,None]*np.column_stack([S,np.ones(self.size)])
        rhs/=storage
        iso=self.cfg.isothermal_after>0 and t>=self.cfg.isothermal_after
        if iso:rhs[:,0]=0
        if not jac:return rhs.ravel()
        bl=np.empty((len(l),2,2));br=np.empty_like(bl)
        bl[:,0,0]=g*kf;br[:,0,0]=-g*kf
        bl[:,0,1]=g*2*(kr/ks)**2*kc[l]*dT;br[:,0,1]=g*2*(kl/ks)**2*kc[r]*dT
        bl[:,1,0]=mtl;br[:,1,0]=mtr;bl[:,1,1]=mcl;br[:,1,1]=mcr
        if iso:bl[:,0,:]=0;br[:,0,:]=0
        vals=[]
        for dest,sgn in ((l,-1.),(r,1.)):
            for deriv in (bl,br):
                for a in range(2):
                    for b in range(2):vals.append(sgn*deriv[:,a,b]/storage[dest,a])
        vals.extend([-rhs[:,0]*Sc/S,-self.bound*self.cfg.hT/storage[:,0]*(0 if iso else 1),-self.bound*self.cfg.hm/self.w])
        return coo_matrix((np.concatenate(vals),(self.jrows,self.jcols)),shape=(2*self.size,2*self.size)).tocsc()
    def profile(self,y):
        field=np.asarray(y).reshape(self.shape+(2,))[:,0,:]
        return PchipInterpolator(self.r**2,field,axis=0)(np.linspace(0,self.R,21)**2)


def derivative_dense(poly,t):
    if hasattr(poly,'D'):
        x=(t-poly.t_shift)/poly.denom;prod=1.;der=0.;ans=np.zeros(poly.D.shape[1])
        for j in range(poly.order):
            der=der*x[j]+prod/poly.denom[j];prod*=x[j];ans+=poly.D[j+1]*der
        return ans
    # Radau dense polynomial y0 + Q*[x,x^2,...]
    x=(t-poly.t_old)/poly.h
    return poly.Q@((np.arange(poly.order+1)+1)*x**np.arange(poly.order+1))/poly.h


def solve(cfg:Config,path:Path,input_path:Path):
    if path.with_suffix('.json').exists(): raise FileExistsError(f'Refusing overwrite: {path}')
    data=load_input(input_path);env=Environment(data,cfg.scenario,cfg.interpolation);op=Operator(cfg)
    y=np.tile([28.,2.55],op.size);times=[0.];samples=[op.profile(y)];snapt=[0.];snaps=[y.copy()]
    fulls=[y.copy()] if cfg.full_samples else []
    event_times=[];event_fields=[]
    diagnostics=[];balances=[];stats={'config':asdict(cfg),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'future':env.future.tolist(),'nodes':op.size,'min_dr_m':float(np.diff(op.r).min()),'max_dr_m':float(np.diff(op.r).max()),'nfev':0,'njev':0,'nlu':0,'steps':0}
    cuts=np.unique(np.r_[data[:,0],np.arange(21600.,cfg.end_s,21600.),cfg.end_s,([cfg.isothermal_after] if cfg.isothermal_after else [])]);cuts=cuts[cuts<=cfg.end_s]
    tic=time.perf_counter();quad=np.polynomial.legendre.leggauss(4);qint=jint=stored=0.;event=None;projections=[]
    minC=2.55;maxC=2.55;minT=28.;maxT=28.;max_up=0.;max_radial=0.;max_axial=0.;max_argmax=0
    root_slope=0.;surface_event=None;mean_event=None
    nextout=60.
    for a,b in zip(cuts[:-1],cuts[1:]):
        if cfg.isothermal_after and a==cfg.isothermal_after:
            dy=np.max(abs(y[::2]-env.future[0]));projections.append({'t_s':float(a),'max_temperature_projection_K':float(dy)})
            y[::2]=env.future[0]
        e=env.segment(a,b)
        fun=lambda t,v:op.rhs(t,v,e)
        def evt(t,v):return np.max(v[1::2])-.15
        evt.direction=-1;evt.terminal=False
        sol=solve_ivp(fun,(a,b),y,method=cfg.method,jac=lambda t,v:op.rhs(t,v,e,True),rtol=cfg.rtol,atol=np.tile([cfg.atol_T,cfg.atol_C],op.size),max_step=cfg.early_step if a<14400 else cfg.max_step,dense_output=True,events=evt)
        if not sol.success:raise RuntimeError(f'Failed on {(a,b)}: {sol.message}')
        if not np.isfinite(sol.y).all() or np.any(sol.y[1::2]<=0):raise FloatingPointError('Invalid accepted field')
        stats['steps']+=len(sol.t)-1
        for k in ('nfev','njev','nlu'):stats[k]+=int(getattr(sol,k))
        if len(sol.t_events[0]):
            root=float(sol.t_events[0][0]);delta=.01
            # Bracket width is reporting evidence, NOT a bound for discretization error.
            tm=root-delta;tp=root+delta
            yr=sol.sol(root);ym=sol.sol(tm);yp=sol.sol(tp)
            slope=float(fun(root,yr)[1::2][np.argmax(yr[1::2])])
            event={'critical_s':root,'critical_h':root/3600,'critical_max':float(np.max(yr[1::2])),'t_minus_s':tm,'max_minus':float(np.max(ym[1::2])),'t_plus_s':tp,'max_plus':float(np.max(yp[1::2])),'critical_slope_kgkg_s':slope,'bracket_s':2*delta,'argmax_flat_index':int(np.argmax(yr[1::2]))}
            actual_end=cfg.execution_time_s if cfg.execution_time_s else tp
            if not (actual_end>root and actual_end<=b):raise ValueError('Execution time must be after the root and in its resolved segment')
            yexec=sol.sol(actual_end)
            event['execution_s']=float(actual_end);event['execution_h']=float(actual_end/3600)
            event['execution_max']=float(np.max(yexec[1::2]))
            event['strict_execution_pass']=bool(np.max(yexec[1::2])<.15)
            if not event['strict_execution_pass']:raise AssertionError('Execution state not strictly qualified')
            event_times=[tm,root,tp,actual_end];event_fields=[ym.copy(),yr.copy(),yp.copy(),yexec.copy()]
        else:actual_end=b
        qtimes=np.arange(nextout,actual_end+1e-8,60.)
        if len(qtimes):
            vals=sol.sol(qtimes)
            for i,t in enumerate(qtimes):
                samples.append(op.profile(vals[:,i]));times.append(float(t))
                if cfg.full_samples:fulls.append(vals[:,i].copy())
            nextout=float(qtimes[-1]+60)
        # Full mesh checks at each accepted time and each mid-step, not just 21 exported radii.
        ct=np.unique(np.r_[sol.t[sol.t<=actual_end],(sol.t[:-1]+sol.t[1:])/2,actual_end]);ct=ct[ct<=actual_end]
        for t in ct:
            u=sol.sol(t).reshape(op.shape+(2,));cc=u[:,:,1];tt=u[:,:,0]
            minC=min(minC,float(cc.min()));maxC=max(maxC,float(cc.max()));minT=min(minT,float(tt.min()));maxT=max(maxT,float(tt.max()))
            max_radial=max(max_radial,float(np.diff(cc,axis=0).max()));max_argmax=max(max_argmax,int(np.argmax(cc)))
            if cfg.nz:max_axial=max(max_axial,float(np.diff(cc,axis=1).max()))
        for i,t in enumerate(sol.t):
            if t>actual_end:break
            u=sol.y[:,i].reshape(op.size,2);mean=float(op.w@u[:,1]/op.vol)
            if diagnostics:max_up=max(max_up,float(np.max(u[:,1])-diagnostics[-1][1]))
            diagnostics.append([float(t),float(np.max(u[:,1])),float(np.min(u[:,1])),mean,float(u[0,0]),float(u[-1,0]),float(np.sum(op.bound*cfg.hm*(u[:,1]-e(t)[1]))/op.vol)])
        # Surface and weighted-average threshold events are diagnostic only.
        for typ,func in [('surface',lambda v:v.reshape(op.shape+(2,))[-1,0,1]-.15),('mean',lambda v:op.w@v[1::2]/op.vol-.15)]:
            if (surface_event if typ=='surface' else mean_event) is None:
                if func(sol.y[:,0])>=0 and func(sol.y[:,-1])<0:
                    tr=brentq(lambda t:func(sol.sol(t)),a,b,xtol=1e-8)
                    if typ=='surface':surface_event=tr
                    else:mean_event=tr
        if cfg.quadrature:
            for lo,hi,poly in zip(sol.t[:-1],sol.t[1:],sol.sol.interpolants):
                if lo>=actual_end:break
                hi=min(hi,actual_end)
                for g,wg in zip(*quad):
                    tq=(lo+hi)/2+(hi-lo)/2*g;wt=wg*(hi-lo)/2
                    uu=poly(tq).reshape(op.size,2);du=derivative_dense(poly,tq).reshape(op.size,2)
                    S=props(uu[:,0],uu[:,1])[0];et,ec=e(tq)
                    jrate=np.sum(op.bound*cfg.hm*(uu[:,1]-ec))/op.vol;qrate=np.sum(op.bound*cfg.hT*(uu[:,0]-et))/op.vol
                    jint+=wt*jrate;qint+=wt*qrate;stored+=wt*(op.w@(S*du[:,0]))/op.vol
            endy=sol.sol(actual_end).reshape(op.size,2)
            balances.append([actual_end,float(op.w@endy[:,1]/op.vol-2.55+jint),stored+qint,jint,qint,stored])
        if b in (1800,3600,5400,7200,9000,10800,14400,21600) or b%21600==0:
            if b<=actual_end:snapt.append(float(b));snaps.append(sol.y[:,-1].copy())
        y=sol.sol(actual_end)
        if event:
            times.append(actual_end);samples.append(op.profile(y));snapt.append(actual_end);snaps.append(y.copy())
            if cfg.full_samples:fulls.append(y.copy())
            break
        if b%21600==0 or b==14400:print(f'n={cfg.n} nz={cfg.nz} {cfg.flux} t={b/3600:g}h M={y[1::2].max():.8f} elapsed={time.perf_counter()-tic:.1f}s',flush=True)
        # Free the previous segment's dense history before allocating the next one.
        sol=None
    stats.update({'elapsed_s':time.perf_counter()-tic,'event':event,'global_C_min':minC,'global_C_max':maxC,'global_T_min':minT,'global_T_max':maxT,'max_radial_increase':max_radial,'max_axial_increase':max_axial,'max_argmax_flat_index':max_argmax,'max_M_increase_accepted':max_up,'surface_crossing_s':surface_event,'mean_crossing_s':mean_event,'projections':projections})
    if balances:
        bb=np.array(balances);stats['max_water_balance_kgkg']=float(abs(bb[:,1]).max());stats['max_effective_heat_balance_J_m3']=float(abs(bb[:,2]).max());stats['effective_heat_relative']=float(abs(bb[:,2]).max()/max(1.,abs(qint)))
    result={'time_s':np.array(times),'radius_cm':np.linspace(0,2,21),'sample_TC':np.array(samples),'r_m':op.r,'z_m':op.z,'weights':op.w.reshape(op.shape),'snapshot_time_s':np.array(snapt),'snapshots':np.array(snaps).reshape((-1,)+op.shape+(2,)),'diagnostics':np.array(diagnostics),'balance':np.array(balances)}
    result['event_time_s']=np.array(event_times)
    result['event_fields_TC']=np.array(event_fields).reshape((-1,)+op.shape+(2,))
    stats['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    stats['input_sha256']=hashlib.sha256(input_path.read_bytes()).hexdigest()
    if cfg.full_samples:result['full_TC']=np.array(fulls).reshape((-1,)+op.shape+(2,))
    path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path.with_suffix('.npz'),**result)
    path.with_suffix('.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(stats,ensure_ascii=False),flush=True)
    return stats

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path(__file__).resolve().parents[1]/'inputs/attachment1.xlsx');p.add_argument('--out',type=Path,required=True);p.add_argument('--config',type=Path)
    p.add_argument('--n',type=int,default=256);p.add_argument('--nz',type=int,default=0);p.add_argument('--flux',default='integral');p.add_argument('--scenario',default='mean');p.add_argument('--rtol',type=float,default=2e-9);p.add_argument('--max-step',type=float,default=240.);p.add_argument('--quadrature',action='store_true');p.add_argument('--full-samples',action='store_true')
    args=p.parse_args();cfg=Config(**json.loads(args.config.read_text())) if args.config else Config(n=args.n,nz=args.nz,flux=args.flux,scenario=args.scenario,rtol=args.rtol,max_step=args.max_step,quadrature=args.quadrature,full_samples=args.full_samples)
    solve(cfg,args.out,args.input)
if __name__=='__main__':main()
