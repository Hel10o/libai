#!/usr/bin/env python3
"""Q1 cylindrical effective transport model; no hidden notebook state.

Units: SI internally, temperature in degC (only temperature differences).
Run from the package root:
  python source/q1_solver.py --input input/附件1.xlsx --out output --n 10240
The validation driver calls the same spatial operator, but independent Bessel
and Chebyshev references are implemented in q1_validate.py and q1_refine.py.
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
from scipy.sparse import diags, csr_matrix

@dataclass(frozen=True)
class Parameters:
    R: float = 0.02
    L: float = 0.25
    rho: float = 820.0
    cp: float = 2600.0
    k: float = 0.36
    hT: float = 25.0
    hm: float = 8e-7
    D0: float = 7e-9
    a: float = 0.89
    T0: float = 28.0
    C0: float = 2.55
    @property
    def alpha(self): return self.k/(self.rho*self.cp)
P=Parameters()
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

def read_xlsx(path: str|Path):
    """Read numeric input via standard-library OOXML, without changing the file.
    Supports shared strings / inline strings / numeric cells; rejects formulas
    without cached numeric values and requires the expected three input columns.
    This independent reader is cross-checked against artifact_tool in this run.
    """
    path=Path(path)
    if not path.is_file(): raise FileNotFoundError(path)
    with zipfile.ZipFile(path) as z:
        bad=z.testzip()
        if bad: raise ValueError(f'ZIP CRC failed: {bad}')
        ss=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS):
                ss.append(''.join(t.text or '' for t in si.iter('{'+NS['s']+'}t')))
        root=ET.fromstring(z.read('xl/workbook.xml'))
        sheets=root.find('s:sheets',NS)
        names=[s.attrib['name'] for s in sheets]
        rid=sheets[0].attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        target=next(e.attrib['Target'] for e in rels if e.attrib['Id']==rid)
        target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
        rows=[]
        for row in ET.fromstring(z.read(target)).findall('s:sheetData/s:row',NS):
            vals={}
            for c in row.findall('s:c',NS):
                col=''.join(ch for ch in c.attrib['r'] if ch.isalpha())
                v=c.find('s:v',NS); typ=c.attrib.get('t')
                if typ=='s': val=ss[int(v.text)] if v is not None else None
                elif typ=='inlineStr': val=''.join(e.text or '' for e in c.findall('.//s:t',NS))
                elif v is None: val=None
                else:
                    try: val=float(v.text)
                    except (ValueError,TypeError): val=v.text
                vals[col]=val
            if vals: rows.append([vals.get(c) for c in ['A','B','C']])
    if rows[0]!=['时间','温度','水分浓度']:
        raise ValueError(f'Unexpected headings: {rows[0]}')
    data=np.asarray(rows[1:],dtype=float)
    if data.ndim!=2 or data.shape[1]!=3 or not np.isfinite(data).all():
        raise ValueError('Missing/nonfinite/non-numeric environment data')
    dt=np.diff(data[:,0])
    if not np.all(dt>0): raise ValueError('Input times must be strictly increasing; duplicate or reversed rows')
    if data[0,0]>0 or data[-1,0]<1800: raise ValueError('Input does not cover [0,1800] s')
    if np.any(data[:,2]<0): raise ValueError('Negative environmental concentration')
    raw=path.read_bytes()
    audit={'input_name':path.name,'sha256':hashlib.sha256(raw).hexdigest(),
           'git_blob_sha1':hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest(),
           'sheet_names':names,'headers':rows[0],'data_rows':len(data),
           'time_min_s':float(data[0,0]),'time_max_s':float(data[-1,0]),
           'time_steps_s':np.unique(dt).tolist(),'missing':0,'duplicate_times':0,
           'temperature_minmax_degC':[float(data[:,1].min()),float(data[:,1].max())],
           'environment_C_minmax':[float(data[:,2].min()),float(data[:,2].max())],
           'temperature_down_steps':int(np.sum(np.diff(data[:,1])<0)),
           'environment_C_down_steps':int(np.sum(np.diff(data[:,2])<0)),
           'interpretation':'External forcing only, not internal response observations',
           'units_source':'Problem appendix 1; xlsx headings contain no unit strings'}
    return data,audit

class Environment:
    def __init__(self,data,kind='linear'):
        self.data=np.asarray(data,dtype=float); self.kind=kind
        self.knots=self.data[:,0]
        if kind not in ('linear','pchip'): raise ValueError(kind)
        self.pchip=PchipInterpolator(self.knots,self.data[:,1:],axis=0,extrapolate=False) if kind=='pchip' else None
    def __call__(self,t):
        t=np.asarray(t)
        if np.any(t<self.knots[0]-1e-10) or np.any(t>self.knots[-1]+1e-10):
            raise ValueError('Environmental extrapolation is prohibited in Q1')
        if self.kind=='pchip': return self.pchip(t)
        return np.stack([np.interp(t,self.knots,self.data[:,i]) for i in (1,2)],axis=-1)

class ConstantEnvironment:
    def __init__(self,T=28.,C=2.55,end=1800.):
        self.T=T; self.C=C; self.knots=np.array([0.,end]); self.kind='constant'
    def __call__(self,t):
        t=np.asarray(t)
        return np.stack([np.full_like(t,self.T,dtype=float),np.full_like(t,self.C,dtype=float)],axis=-1)

class RadialFV:
    """Vertex-centred finite volumes with nodes at r=0 and r=R.
    w=integral r dr over each dual cell; actual volume=2*pi*L*w.
    Boundary nodes have their proper half-cell storage, not ghost/centre output.
    The final augmented state integrates outward flux divided by total volume.
    """
    def __init__(self,n,field,p=P,transfer_scale=1.,constant_D=None):
        if n<4: raise ValueError('Need n>=4 radial intervals')
        self.p=p; self.n=n; self.field=field; self.constant_D=constant_D
        self.r=np.linspace(0,p.R,n+1); self.dr=p.R/n
        self.faces=np.r_[0.,(self.r[:-1]+self.r[1:])/2,p.R]
        self.w=np.diff(self.faces**2)/2; self.vol=self.w.sum()
        self.s=self.faces[1:-1]/self.dr
        if field=='T': self.beta=p.hT/(p.rho*p.cp)*transfer_scale; self.col=0
        elif field=='C': self.beta=p.hm*transfer_scale; self.col=1
        else: raise ValueError(field)
    def diff(self,u):
        if self.field=='T': return np.full_like(u,self.p.alpha),np.zeros_like(u)
        if self.constant_D is not None: return np.full_like(u,self.constant_D),np.zeros_like(u)
        # No clipping of the solution: a nonpositive nonlinear iterate is a failure.
        if not np.all(np.isfinite(u)) or np.any(u<=0):
            raise FloatingPointError('Nonpositive/nonfinite C in nonlinear solver')
        D=self.p.D0*np.exp(-self.p.a/u)
        return D,D*self.p.a/u**2
    def rhs(self,t,y,env):
        u=y[:-1]; D,_=self.diff(u)
        df=2*D[:-1]*D[1:]/(D[:-1]+D[1:])
        g=self.s*df*(u[:-1]-u[1:]) # positive outward
        f=np.zeros_like(y)
        f[:-2]-=g/self.w[:-1];f[1:-1]+=g/self.w[1:]
        go=self.p.R*self.beta*(u[-1]-env(t)[self.col])
        f[-2]-=go/self.w[-1];f[-1]=go/self.vol
        return f
    def jac(self,t,y,env):
        u=y[:-1];D,Dp=self.diff(u);sumD=D[:-1]+D[1:]
        df=2*D[:-1]*D[1:]/sumD
        dleft=2*(D[1:]/sumD)**2*Dp[:-1]
        dright=2*(D[:-1]/sumD)**2*Dp[1:]
        jump=u[:-1]-u[1:]
        a=self.s*(df+jump*dleft);b=self.s*(-df+jump*dright)
        m=len(u);diag=np.zeros(m+1)
        diag[:m-1]-=a/self.w[:-1];diag[1:m]+=b/self.w[1:]
        outder=self.p.R*self.beta
        diag[m-1]-=outder/self.w[-1]
        lower=np.r_[a/self.w[1:],outder/self.vol]
        upper=np.r_[-b/self.w[:-1],0.]
        return diags([lower,diag,upper],[-1,0,1],format='csc')

def solve_fv(env,field,n=640,p=P,rtol=2e-10,atol=2e-12,max_step=15.,
             end=1800.,times=None,initial=None,transfer_scale=1.,constant_D=None,
             quadrature=False,verbose=False):
    op=RadialFV(n,field,p,transfer_scale,constant_D)
    times=np.arange(0.,end+0.5,1.) if times is None else np.asarray(times,dtype=float)
    if times[0]!=0 or times[-1]!=end: raise ValueError('Requested times must include 0 and end')
    start=p.T0 if field=='T' else p.C0
    if initial is None: u0=np.full(n+1,start)
    elif callable(initial): u0=np.asarray(initial(op.r),dtype=float)
    else: u0=np.broadcast_to(initial,(n+1,)).astype(float).copy()
    y=np.r_[u0,0.]
    rout=np.linspace(0,p.R,21)
    # All production grids have intervals divisible by 20, so no spatial interpolation.
    if n%20: raise ValueError('n must be divisible by 20 for exact output-node selection')
    oi=np.arange(21)*(n//20)
    out=np.empty((len(times),21));out[0]=u0[oi]
    avg=np.empty(len(times));avg[0]=op.w@u0/op.vol
    acc=np.zeros(len(times));int_error=[]
    snapshots={0.:u0.copy()}; selected={100.,300.,600.,900.,1200.,1500.,1800.}
    cuts=np.unique(np.r_[0.,env.knots[(env.knots>0)&(env.knots<end)],end])
    stats={'field':field,'n_intervals':n,'n_nodes':n+1,'dr_m':op.dr,'rtol':rtol,'atol':atol,
           'max_step_s':max_step,'method':'scipy.solve_ivp BDF, analytic sparse Jacobian',
           'segments':len(cuts)-1,'nfev':0,'njev':0,'nlu':0,'accepted_steps':0,
           'nonlinear_failure_policy':'raise; never silently accept or clip',
           'nonlinear_iteration':'SciPy BDF Newton iteration; maximum 4 iterations; unsuccessful convergence reduces step size',
           'first_step_policy':'C only: min(segment length, 1e-3 s, 0.05*dr^2/D0); avoids negative initial-step probe',
           'newton_tolerance':max(10*np.finfo(float).eps/rtol,min(0.03,np.sqrt(rtol))),
           'environment_interpolation':env.kind}
    tic=time.perf_counter(); cumulative_quad=0.
    gx,gw=np.polynomial.legendre.leggauss(4)
    for left,right in zip(cuts[:-1],cuts[1:]):
        sol=solve_ivp(lambda t,y:op.rhs(t,y,env),(left,right),y,
                      method='BDF',jac=lambda t,y:op.jac(t,y,env),rtol=rtol,atol=atol,
                      max_step=max_step,dense_output=True,
                      first_step=min(right-left,1e-3,0.05*op.dr**2/p.D0) if field=='C' else None)
        if not sol.success: raise RuntimeError(f'{field}, n={n}, [{left},{right}]: {sol.message}')
        for key in ('nfev','njev','nlu'): stats[key]+=int(getattr(sol,key))
        stats['accepted_steps']+=len(sol.t)-1
        mask=(times>left)&(times<=right+1e-9)
        ts=times[mask]; ys=sol.sol(ts)
        out[mask]=ys[oi].T
        avg[mask]=(op.w@ys[:-1])/op.vol;acc[mask]=ys[-1]
        for t in selected:
            if left<t<=right: snapshots[t]=sol.sol(t)[:-1].copy()
        if quadrature:
            for k in range(0,len(sol.t)-1,32):
                ta=sol.t[k:min(k+32,len(sol.t)-1)];tb=sol.t[k+1:min(k+33,len(sol.t))]
                tq=((ta+tb)[:,None]/2+(tb-ta)[:,None]/2*gx).ravel()
                vals=sol.sol(tq)[-2]
                flow=op.p.R*op.beta*(vals-env(tq)[:,op.col])/op.vol
                cumulative_quad+=float(np.sum(flow.reshape(-1,4)*gw*(tb-ta)[:,None]/2))
            int_error.append([right,float(op.w@sol.y[:-1,-1]/op.vol-avg[0]+cumulative_quad)])
        y=sol.y[:,-1]
        if verbose and int(right)%300==0: print(f'{field} n={n} t={right:.0f}',flush=True)
    stats['elapsed_seconds']=time.perf_counter()-tic
    stats['conservation_augmented_max_abs']=float(np.max(np.abs(avg-avg[0]+acc)))
    stats['independent_flux_quadrature_max_abs']=float(np.max(np.abs(np.array(int_error)[:,1]))) if int_error else None
    stats['min_value']=float(out.min());stats['max_value']=float(out.max())
    stats['initial_max_error']=float(np.max(np.abs(out[0]-u0[oi])))
    return {'time_s':times,'radius_m':rout,'values':out,'average':avg,'outward_flux_integral_average_units':acc,
            'internal_radius_m':op.r,'internal_weights_rdr':op.w,
            'snapshot_times_s':np.array(sorted(snapshots)),
            'internal_snapshots':np.array([snapshots[t] for t in sorted(snapshots)]),
            'quadrature_errors':np.array(int_error),'stats':stats}

def save_field(path,result):
    np.savez_compressed(path,**{k:v for k,v in result.items() if k!='stats'})
    Path(str(path)+'.json').write_text(json.dumps(result['stats'],indent=2,ensure_ascii=False),encoding='utf-8')

def write_csvs(out,T,C):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for name,res in [('temperature',T),('moisture',C)]:
        np.savetxt(out/f'{name}_unrounded.csv',np.column_stack([res['time_s'],res['values']]),
                   delimiter=',',fmt='%.17g',header='time_s,'+','.join(f'r_{r*100:.1f}_cm' for r in res['radius_m']),comments='')
        ii=[np.flatnonzero(res['time_s']==t)[0] for t in (100,300,600,900,1200,1500,1800)]
        np.savetxt(out/f'table_{name}.csv',np.column_stack([res['time_s'][ii],res['values'][ii][:,::5]]),
                   delimiter=',',fmt=['%.0f']+['%.4f']*5,header='time_s,0_cm,0.5_cm,1_cm,1.5_cm,2_cm',comments='')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=Path(__file__).resolve().parents[1]/'input'/'附件1.xlsx')
    ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'output')
    ap.add_argument('--n',type=int,default=10240)
    ap.add_argument('--rtol',type=float,default=2e-12);ap.add_argument('--atol',type=float,default=2e-14)
    ap.add_argument('--max-step',type=float,default=3.)
    ap.add_argument('--interp',choices=['linear','pchip'],default='linear')
    ap.add_argument('--quadrature',action='store_true')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    data,audit=read_xlsx(args.input);env=Environment(data,args.interp)
    (args.out/'input_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    results={}
    for field in ('T','C'):
        results[field]=solve_fv(env,field,n=args.n,rtol=args.rtol,atol=args.atol,
                                max_step=args.max_step,quadrature=args.quadrature,verbose=True)
        save_field(args.out/f'{field}_n{args.n}.npz',results[field])
        print(json.dumps(results[field]['stats'],ensure_ascii=False),flush=True)
    T,C=results['T'],results['C']
    np.savez_compressed(args.out/'q1_unrounded.npz',time_s=T['time_s'],radius_cm=T['radius_m']*100,
                        temperature_degC=T['values'],moisture_dry_basis=C['values'],
                        average_temperature_degC=T['average'],average_C=C['average'],
                        environment=env(T['time_s']))
    write_csvs(args.out,T,C)
    run={'parameters':asdict(P),'settings':vars(args)|{'input':str(args.input),'out':str(args.out)},
         'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
         'T':T['stats'],'C':C['stats']}
    (args.out/'run.json').write_text(json.dumps(run,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Tables:')
    for field in ('T','C'): print(field,results[field]['values'][[100,300,600,900,1200,1500,1800]][:,::5])

if __name__=='__main__': main()
