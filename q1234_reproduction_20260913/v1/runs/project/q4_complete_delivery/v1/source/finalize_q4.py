"""Select a strict execution point using computed numerical discrepancies.
All reported cells, tables and trajectories originate from the same main run.
The budget is empirical numerical evidence, not a rigorous error bound.
"""
import argparse,json,math,csv
from pathlib import Path
import numpy as np
from q4_common import ROOT

def finalize(base):
    base=Path(base);out=base/'output';v=base/'validation';g=v/'geometry'
    if (out/'result4_data.json').exists():raise FileExistsError('Refuse overwrite final derived results')
    read=lambda p:json.loads(Path(p).read_text())
    main=read(out/'main.json');root=main['event']['critical_s']
    a=read(v/'spectral80.json')['event']['critical_s']
    b=read(v/'time120_bdf.json')['event']['critical_s']
    f320=read(g/'q4_320x0_integral_rtol2e-10.json')['event_root_s']
    f640=read(g/'q4_640x0_integral_rtol2e-10.json')['event_root_s']
    floose=read(g/'q4_320x0_integral.json')['event_root_s']
    fex=(4*f640-f320)/3
    discrepancies={'spectral80_to_main_s':abs(a-root),'same_n_time_method_change_s':abs(b-root),
        'independent_fv_extrapolation_to_main_s':abs(fex-root),
        'independent_fv_temporal_change_s':abs(f320-floose),
        'independent_fv_last_spatial_correction_s':abs(f640-fex)}
    budget=math.ceil(2000*max(discrepancies.values()))/1000
    execute=.36*math.ceil((root+budget)/.36)
    z=dict(np.load(out/'main.npz'));candidates=z['event_time_s'];idx=int(np.argmin(abs(candidates-execute)))
    if abs(candidates[idx]-execute)>1e-7:raise ValueError('Execution point not actually solved')
    full=z['event_full_TC'][idx]
    # Current implementation records exact candidates around the event. A shifted
    # final endpoint must first be regenerated from that complete state.
    if abs(z['time_s'][-1]-execute)>1e-7:raise ValueError('Regenerate final row from event full state before export')
    if float(full[:,1].max())>=.15:raise AssertionError('Strict C<.15 failed')
    before=z['event_full_TC'][np.argmin(abs(candidates-(root-1)))]
    after=z['event_full_TC'][np.argmin(abs(candidates-(root+1)))]
    final={'critical_s':root,'critical_h':root/3600,'critical_display_h':round(root/3600,4),
       'execution_s':execute,'execution_h':execute/3600,'execution_max_C':float(full[:,1].max()),
       'actual_margin_s':execute-root,'empirical_budget_s':budget,'budget_rule':'ceil(2*max(recorded discrepancies), 0.001 s)',
       'budget_is_rigorous_bound':False,'discrepancies':discrepancies,'fv_extrapolated_root_s':fex,
       'one_second_before_max_C':float(before[:,1].max()),'one_second_after_max_C':float(after[:,1].max()),
       'radius_m':float(z['radius_m'][-1]),'length_m':.25,'domain':'0<=r<=R(t); -0.125<=z<=0.125 m',
       'source':'output/main.npz; official 1D endpoint, matched 2D applicability and whole-domain evidence are reported separately'}
    (out/'end_event.json').write_text(json.dumps(final,ensure_ascii=False,indent=2))
    keep=z['time_s']>0;Cf=z['sample_TC'][keep,:,1]
    doc={'time_s':z['time_s'][keep].tolist(),'radius_m':z['radius_m'][keep].tolist(),
       'C_fixed':[[None if not np.isfinite(x) else float(x) for x in row] for row in Cf],
       'C_surface':z['surface_TC'][keep,1].tolist()}
    (out/'result4_data.json').write_text(json.dumps(doc,ensure_ascii=False,allow_nan=False))
    with (out/'trajectory_60s_unrounded.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','radius_m','max_C','max_radius_m','mean_C','T_surface_C','C_surface']+
          [f'T_r{j/10:.1f}cm_C' for j in range(21)]+[f'C_r{j/10:.1f}cm' for j in range(21)])
        for i,t in enumerate(z['time_s']):
            vals=[t,z['radius_m'][i],z['max_C'][i],z['max_radius_m'][i],z['mean_C'][i],*z['surface_TC'][i],*z['sample_TC'][i,:,0],*z['sample_TC'][i,:,1]]
            w.writerow(['' if not np.isfinite(x) else format(float(x),'.17g') for x in vals])
    rows=[]
    for t in np.r_[np.arange(21600,execute,21600),execute]:
        i=int(np.argmin(abs(z['time_s']-t)))
        if abs(z['time_s'][i]-t)>1e-7:raise ValueError('Missing actual table time')
        rows.append([t/3600,*z['sample_TC'][i,::5,1],z['surface_TC'][i,1],z['radius_m'][i]*100])
    with (out/'table6_unrounded.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_h','r0cm','r0.5cm','r1cm','r1.5cm','r2cm','surface','radius_cm'])
        for row in rows:w.writerow(['' if not np.isfinite(x) else format(float(x),'.17g') for x in row])
    md='|时间/h|0 cm|0.5 cm|1 cm|1.5 cm|2 cm|药材表面|半径/cm|\n|---:|---:|---:|---:|---:|---:|---:|---:|\n'
    for row in rows:md+='|'+'|'.join('—' if not np.isfinite(x) else f'{x:.4f}' for x in row)+'|\n'
    md+='\n“—”表示固定位置已在当前材料域外。末行轴心显示0.1500是四位舍入；未舍入最大值见end_event.json，严格小于0.15。\n'
    (out/'表6.md').write_text(md)
    print(json.dumps(final,ensure_ascii=False,indent=2));print(md)
    return final
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=ROOT);a=p.parse_args();finalize(a.base)
