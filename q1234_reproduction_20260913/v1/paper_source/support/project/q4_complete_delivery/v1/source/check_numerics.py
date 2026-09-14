"""Read-only comparisons of real saved runs. No output-state edits."""
import argparse,json,hashlib
from pathlib import Path
from decimal import Decimal,ROUND_HALF_UP
import numpy as np
from q4_common import ROOT
from q4_spectral import interp,primitive,inverse_primitive

def check(base=ROOT):
    base=Path(base);v=base/'validation';z=np.load(base/'output/main.npz');m=json.loads((base/'output/main.json').read_text())
    result={'main':'output/main.npz','comparisons':{},'main_mass_balance_max_mean_C':m['mass_balance_max_abs_in_mean_C'],
       'main_max_boundary_residual':m['max_boundary_residual'],'main_total_loss_per_unit_dry_mass':float(2.55-z['mean_C'][-1]),
       'main_end_mean_C':float(z['mean_C'][-1]),'domain_masks':{}}
    for name in ['spectral80','spectral160','time120_bdf']:
        p=v/(name+'.npz')
        if not p.exists():result['comparisons'][name]={'status':'not complete; excluded'};continue
        q=np.load(p);d=json.loads(p.with_suffix('.json').read_text());common,ia,ib=np.intersect1d(z['time_s'],q['time_s'],return_indices=True)
        err=z['sample_TC'][ia]-q['sample_TC'][ib];finite=np.isfinite(err)
        rnd=lambda x:Decimal(str(float(x))).quantize(Decimal('.0001'),rounding=ROUND_HALF_UP)
        mismatch=sum(rnd(a)!=rnd(b) for a,b in zip(z['sample_TC'][ia][finite],q['sample_TC'][ib][finite]))
        result['comparisons'][name]={'status':'actual saved run read','common_times':len(common),
            'max_abs_T_C':float(np.nanmax(abs(err[:,:,0]))),'max_abs_C':float(np.nanmax(abs(err[:,:,1]))),
            'four_decimal_mismatch_count_TC':int(mismatch),'root_difference_s':d['event']['critical_s']-m['event']['critical_s']}
    g=v/'geometry';p=g/'q4_640x0_integral_rtol2e-11.npz'
    if p.exists():
        f=np.load(p);common,ia,ib=np.intersect1d(z['time_s'],f['time_s'],return_indices=True)
        n=len(z['x'])-1;w=(-1.)**np.arange(n+1);w[[0,-1]]*=.5
        I=interp(z['x'],w,np.linspace(0,1,21)**2);u=z['full_TC'][ia]
        mt=np.einsum('ij,tj->ti',I,u[:,:,0]);mc=inverse_primitive(np.einsum('ij,tj->ti',I,primitive(u[:,:,1])),np.einsum('ij,tj->ti',I,u[:,:,1]))
        result['independent_fv']={'path':str(p.relative_to(base)),'common_times':len(common),'coordinate':'matched material xi=0:.05:1, not fixed physical r',
            'max_abs_T':float(np.max(abs(mt-f['mid'][ib,:,0]))),'max_abs_C':float(np.max(abs(mc-f['mid'][ib,:,1]))),
            'last_pair_extrapolation_is_estimate':True}
    mask=z['fixed_radius_m'][None,:]>z['radius_m'][:,None]+1e-14
    result['domain_masks']={'outside_count':int(mask.sum()),'outside_are_nan':bool(np.isnan(z['sample_TC'][:,:,1][mask]).all()),
        'inside_are_finite':bool(np.isfinite(z['sample_TC'][:,:,1][~mask]).all())}
    eq=np.isclose(z['fixed_radius_m'][None,:],z['radius_m'][:,None],atol=1e-14,rtol=0);ii,jj=np.where(eq)
    result['surface_fixed_coincidence']={'count':len(ii),'max_abs_C_difference':float(np.max(abs(z['sample_TC'][ii,jj,1]-z['surface_TC'][ii,1])))}
    result['polynomial_maximum']={'max_C_minus_axis_C':float(np.max(abs(z['max_C']-z['full_TC'][:,0,1]))),
       'root_extrema':m['event']['root_extrema'],'early_argmax_note':'nearly uniform profiles have floating-point nonunique locations'}
    (v/'numeric_comparison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False,indent=2))
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=ROOT);a=p.parse_args();check(a.base)
