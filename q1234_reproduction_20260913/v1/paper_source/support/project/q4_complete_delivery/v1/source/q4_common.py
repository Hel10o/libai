"""Inputs and geometry, all in SI. Only reads source workbooks."""
from pathlib import Path
import hashlib, json
import numpy as np
from openpyxl import load_workbook
from scipy.interpolate import PchipInterpolator

ROOT=Path(__file__).resolve().parents[1]

def read_numeric(path, columns):
    wb=load_workbook(path,read_only=True,data_only=True)
    ws=wb.worksheets[0];rows=list(ws.values);wb.close()
    data=np.asarray([r[:columns] for r in rows[1:] if r[0] is not None],float)
    if not np.isfinite(data).all() or np.any(np.diff(data[:,0])<=0):
        raise ValueError('Nonfinite input, duplicate or unordered time')
    return rows[0],data

def load_input(path):return read_numeric(path,3)[1]

class Environment:
    def __init__(self,data):
        self.data=np.asarray(data);self.future=self.data[self.data[:,0]>=10800,1:].mean(axis=0)
    def __call__(self,t):
        if t>14400:return self.future
        return np.array([np.interp(t,self.data[:,0],self.data[:,k]) for k in (1,2)])
    def segment(self,a,b):return (lambda t:self.future) if a>=14400 else self

class Radius:
    def __init__(self,data,kind='linear'):
        self.data=np.asarray(data).copy();self.data[:,1]*=.01;self.kind=kind
        if np.any(self.data[:,1]<=0) or np.any(np.diff(self.data[:,1])>0):raise ValueError('Invalid measured shrinkage')
        self.pchip=PchipInterpolator(self.data[:,0],self.data[:,1],extrapolate=False)
    def __call__(self,t):
        if self.kind=='fixed':return .02
        if self.kind=='pchip' and 0<=t<=self.data[-1,0]:return float(self.pchip(t))
        return float(np.interp(t,self.data[:,0],self.data[:,1]))
    def derivative(self,t):
        if self.kind=='fixed' or t>=self.data[-1,0]:return 0.
        if self.kind=='pchip':return float(self.pchip.derivative()(t))
        i=np.clip(np.searchsorted(self.data[:,0],t,side='right')-1,0,len(self.data)-2)
        return float(np.diff(self.data[i:i+2,1])[0]/np.diff(self.data[i:i+2,0])[0])

def input_audit(out):
    result={}
    for name,n in [('attachment1.xlsx',3),('attachment2.xlsx',2)]:
        p=ROOT/'inputs'/name;head,a=read_numeric(p,n)
        result[name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'headers':list(head),
          'rows':len(a),'first':a[0].tolist(),'last':a[-1].tolist(),'missing':int(np.isnan(a).sum()),
          'time_step_unique_s':np.unique(np.diff(a[:,0])).tolist(),'min':a.min(axis=0).tolist(),'max':a.max(axis=0).tolist()}
    _,a=read_numeric(ROOT/'inputs/attachment2.xlsx',2)
    result['radius']={'units_input':'cm','units_internal':'m','monotone_nonincreasing':bool(np.all(np.diff(a[:,1])<=0)),
       'plateau_intervals':int(np.sum(np.diff(a[:,1])==0)),'measurement_increment_cm':.001,
       'primary_interpolation':'piecewise linear, no overshoot; restart at knots',
       'beyond_72h':'hold last observed radius; conditional extension only if needed',
       'length_m':.25,'length_assumption':'constant; homogeneous radial material contraction'}
    result['future_environment']=Environment(load_input(ROOT/'inputs/attachment1.xlsx')).future.tolist()
    Path(out).write_text(json.dumps(result,ensure_ascii=False,indent=2))
    return result
