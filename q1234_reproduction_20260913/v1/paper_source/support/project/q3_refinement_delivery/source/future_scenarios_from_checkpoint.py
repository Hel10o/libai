"""High-resolution N=2048 future-boundary comparisons from the verified 4 h full state.

This script reuses only a full internal checkpoint whose provenance is fixed by SHA256.
It does not reuse the 21-point display profile. From 4 h onward each requested
future environment is integrated with the same FV operator, tolerances, and all-domain
event definition as the historical main model.
"""
from __future__ import annotations
import argparse, hashlib, json, time, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from q3_solver import Config, Operator, load_input


def future_value(data, scenario):
    if scenario == 'mean':
        return data[data[:,0] >= 10800,1:].mean(axis=0)
    if scenario == 'time_mean':
        tail=data[data[:,0] >= 10800]
        return np.trapezoid(tail[:,1:], tail[:,0], axis=0)/(tail[-1,0]-tail[0,0])
    if scenario == 'nominal':
        return np.array([50.0,0.05])
    if scenario == 'last':
        return data[-1,1:].copy()
    raise ValueError(scenario)


def run(checkpoint, input_file, scenario, out):
    if out.exists():
        raise FileExistsError(out)
    data=load_input(input_file)
    z=np.load(checkpoint)
    if 'state_TC' in z.files:
        if abs(float(z['time_s'])-14400.0)>1e-8: raise ValueError('Checkpoint is not at 4 h')
        y=z['state_TC'].reshape(-1)
    else:
        st=z['snapshot_time_s']
        idx=np.where(np.isclose(st,14400.0))[0]
        if len(idx)!=1: raise ValueError('No unique 4 h checkpoint')
        y=z['snapshots'][idx[0]].reshape(-1)
    cfg=Config(n=2048,rtol=2e-12,atol_T=2e-12,atol_C=2e-14,max_step=60,early_step=15,quadrature=False,full_samples=False,execution_time_s=0)
    op=Operator(cfg)
    if y.size != 2*op.size: raise ValueError('Checkpoint grid mismatch')
    fut=future_value(data,scenario)
    env=lambda t: np.broadcast_to(fut,np.shape(t)+(2,))
    fun=lambda t,v: op.rhs(t,v,env)
    def ev(t,v): return float(np.max(v[1::2])-0.15)
    ev.direction=-1;ev.terminal=False
    tic=time.perf_counter()
    sol=solve_ivp(fun,(14400.0,206940.0),y,method='BDF',jac=lambda t,v:op.rhs(t,v,env,True),rtol=cfg.rtol,atol=np.tile([cfg.atol_T,cfg.atol_C],op.size),max_step=cfg.max_step,dense_output=True,events=ev)
    if not sol.success or not len(sol.t_events[0]): raise RuntimeError(sol.message or 'event not reached')
    root=float(sol.t_events[0][0]); yr=sol.sol(root); rootmax=float(np.max(yr[1::2])); argmax=int(np.argmax(yr[1::2])); execution=0.36*math.ceil((root+0.03)/0.36-1e-13); yexec=sol.sol(execution)
    state_206901=sol.sol(206901.0)
    obj={
      'scenario':scenario,'future_T_C':float(fut[0]),'future_C_kgkg':float(fut[1]),
      'critical_s':root,'critical_h':root/3600,'critical_max':rootmax,'critical_argmax_flat_index':argmax,'execution_s':execution,'execution_h':execution/3600,'execution_max':float(np.max(yexec[1::2])),'strict_execution_pass':bool(np.max(yexec[1::2])<0.15),
      'Cmax_206901':None if state_206901 is None else float(np.max(state_206901[1::2])),
      'strict_at_206901':None if state_206901 is None else bool(np.max(state_206901[1::2])<0.15),
      'elapsed_s':time.perf_counter()-tic,'nfev':sol.nfev,'njev':sol.njev,'nlu':sol.nlu,
      'checkpoint_sha256':hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
      'input_sha256':hashlib.sha256(input_file.read_bytes()).hexdigest(),
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(obj,ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--scenario',choices=['mean','time_mean','nominal','last'],required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.checkpoint,a.input,a.scenario,a.out)
