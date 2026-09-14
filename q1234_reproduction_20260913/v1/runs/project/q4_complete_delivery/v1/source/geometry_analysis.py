"""Summarize completed geometry runs; never treats absent runs as passed."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from geometry_solver import DEST

def summarize():
    rows=[];details={};files=[];event_checks={}
    for f in sorted(DEST.glob('q*x*.json')):
        s=json.loads(f.read_text())
        if not s.get('nz'):continue
        suffix='_integral' if s['flux'].startswith('integral') else ''
        base=f.stem.replace(f"{s['case']}_{s['nr']}x{s['nz']}",f"{s['case']}_{s['nr']}x0",1)
        if not (DEST/(base+'.json')).exists():continue
        one=json.loads((DEST/(base+'.json')).read_text())
        a=np.load(DEST/(base+'.npz'));b=np.load(f.with_suffix('.npz'))
        ts,ia,ib=np.intersect1d(a['time_s'],b['time_s'],return_indices=True)
        mid=b['mid'][ib]-a['mid'][ia];end=b['end'][ib]-a['mid'][ia]
        row={'case':s['case'],'nr':s['nr'],'nz':s['nz'],'flux':suffix.strip('_') or 'harmonic','rtol':s['rtol'],'max_step_s':s['max_step'],
             'event_1d_s':one['event_root_s'],'event_2d_s':s['event_root_s'],
             'end_effect_s':None if one['event_root_s'] is None or s['event_root_s'] is None else s['event_root_s']-one['event_root_s'],
             'max_mid_T_difference_all_samples_K':float(abs(mid[:,:,0]).max()),
             'max_mid_C_difference_all_samples':float(abs(mid[:,:,1]).max()),
             'max_end_T_difference_all_samples_K':float(abs(end[:,:,0]).max()),
             'max_end_C_difference_all_samples':float(abs(end[:,:,1]).max()),
             'paired_sample_count':len(ts),'elapsed_s':s['elapsed_s']}
        masks={'q1_table_range':ts<=1800,'q2_table_range':ts<=10800,'full_common_range':np.ones(len(ts),bool)}
        d={'run':f.stem,'paired_1d':base,'summary':row,'intervals':{},'end_location_r_z_m':s['end_max_location_r_z_m']}
        for label,mask in masks.items():
            d['intervals'][label]={'last_sample_s':float(ts[mask][-1]),'max_mid_difference_T_C':np.max(abs(mid[mask]),axis=(0,1)).tolist(),
                 'max_end_difference_T_C':np.max(abs(end[mask]),axis=(0,1)).tolist()}
        if s['event_root_s'] is not None:
            d['event']={'time_s':b['event_times_s'].tolist(),'Cmax':[float(v[:,:,1].max()) for v in b['event_states']],
                        'root_max_r_z_m':s['end_max_location_r_z_m']}
            states=b['event_states'];c=states[:,:,:,1]
            ec={'all_event_fields_finite':bool(np.isfinite(states).all()),'minimum_C':float(c.min()),
                'largest_positive_radial_C_increment':float(np.max(np.diff(c,axis=2))),
                'largest_positive_axial_C_increment':float(np.max(np.diff(c,axis=1))),
                'max_locations_zi_ri':[list(map(int,np.unravel_index(np.argmax(v),v.shape))) for v in c]}
            assert ec['all_event_fields_finite'] and ec['minimum_C']>0
            assert ec['largest_positive_radial_C_increment']<1e-10 and ec['largest_positive_axial_C_increment']<1e-10
            event_checks[f.stem]=ec
        profiles=[]
        for t in [100,300,600,900,1200,1500,1800,3600,5400,7200,9000,10800,21600,43200,86400,129600,172800]:
            hit=np.flatnonzero(ts==t)
            if len(hit):
                j=int(hit[0]);profiles.append({'time_s':t,'reference_radius_fraction':np.linspace(0,1,21).tolist(),
                      'actual_R_m':float(b['radius_m'][ib[j]]),'mid_2d_T_C':b['mid'][ib[j]].tolist(),
                      'end_2d_T_C':b['end'][ib[j]].tolist(),'radial_1d_T_C':a['mid'][ia[j]].tolist()})
        d['profiles']=profiles
        rows.append(row);details[f.stem]=d
        for p in [f,f.with_suffix('.npz'),DEST/(base+'.json'),DEST/(base+'.npz')]:
            files.append({'path':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (DEST/'geometry_comparison.json').write_text(json.dumps(details,indent=2,ensure_ascii=False)+'\n')
    (DEST/'event_field_shape_checks.json').write_text(json.dumps(event_checks,indent=2)+'\n')
    if rows:
        with (DEST/'geometry_comparison.csv').open('w',newline='') as stream:
            w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (DEST/'geometry_evidence_manifest.json').write_text(json.dumps({'completed_runs':len(rows),'files':list({a['path']:a for a in files}.values())},indent=2)+'\n')
    print(json.dumps(rows,ensure_ascii=False,indent=2))
    return rows,details

if __name__=='__main__':summarize()
