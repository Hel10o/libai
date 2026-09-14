"""Read-only numeric/source audit for the newly solved Q2/Q3 trajectory.
Frozen historical comparison inputs are included in inputs/historical_baseline.
All new evidence goes to --out, which must not already exist.
"""
from pathlib import Path
import argparse,json,csv,hashlib,zipfile,xml.etree.ElementTree as ET
import numpy as np
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def fp(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def rounded(a):
    r=np.round(a,4)
    # Resolve numbers close enough to a tie to expose decimal formatting choices.
    z=a*10000;ii=np.argwhere(abs(z-np.floor(z)-.5)<1e-9)
    for ix in ii:
        ix=tuple(ix);r[ix]=float(format(float(a[ix]),'.4f'))
    return r

def cmp(a,b,time_s):
    d=abs(a-b);rr=rounded(a)!=rounded(b)
    report={'max_abs_TC':d.max(axis=(0,1)).tolist(),'rms_TC':np.sqrt(np.mean(d*d,axis=(0,1))).tolist(),'four_decimal_difference_count_TC':rr.sum(axis=(0,1)).tolist(),'values_per_field':int(np.prod(d.shape[:2]))}
    report['argmax_TC']=[]
    for k in range(2):
        i,j=np.unravel_index(d[:,:,k].argmax(),d.shape[:2]);report['argmax_TC'].append({'t_s':float(time_s[i]),'radius_cm':float(j/10),'new':float(a[i,j,k]),'comparison':float(b[i,j,k])})
    return report,rr

def old_workbook(p):
    out=np.empty((10800,21,2));times=np.empty((10800,2))
    with zipfile.ZipFile(p) as z:
        rel={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheets=ET.fromstring(z.read('xl/workbook.xml')).find(NS+'sheets')
        names=[]
        for k,s in enumerate(sheets):
            names.append(s.attrib['name']);target=rel[s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            for event,el in ET.iterparse(z.open(target),events=('end',)):
                if el.tag!=NS+'row':continue
                i=int(el.attrib['r'])-2
                if i<0:el.clear();continue
                if i>=10800:raise ValueError('Unexpected old Q2 row')
                for cell in el:
                    addr=cell.attrib['r'];col=''.join(filter(str.isalpha,addr));v=cell.find(NS+'v')
                    if v is None:continue
                    if col=='A':times[i,k]=float(v.text)
                    elif len(col)==1 and 'B'<=col<='V':out[i,ord(col)-ord('B'),k]=float(v.text)
                el.clear()
    assert names==['温度','水分浓度'] and np.array_equal(times[:,0],np.arange(1,10801)) and np.array_equal(times[:,0],times[:,1])
    return out

def run(repo,new,reference,out):
    out.mkdir(parents=True,exist_ok=False)
    z=np.load(new);t=z['time_s'];v=z['profile_TC'];ref=np.load(reference);rv=ref['profile_TC'];rt=ref['time_s'];ri=np.searchsorted(t,rt);assert np.array_equal(t[ri],rt)
    old_path=repo/'q2_final_delivery/output/q2_unrounded.npz';o=np.load(old_path);ov=np.stack([o['temperature_degC'],o['moisture_dry_basis']],axis=-1)
    r={'source_sha':'anonymous-source','data_contract':{'shape':list(v.shape),'time_s_first_last':[float(t[0]),float(t[-1])],'integer_seconds_0_through':206906,'extra_final_s':206906.76,'all_integer_seconds_present':bool(np.array_equal(t[:-1],np.arange(206907))),'radius_cm':z['radius_cm'].tolist(),'joint_fields':['temperature_degC','C_kg_water_per_kg_dry_solid'],'finite':bool(np.isfinite(v).all())}}
    r['historical_n160_BDF_reference_at_common_outputs'],rr=cmp(v[ri],rv,rt)
    r['historical_n160_BDF_reference_at_common_outputs']['scope']='Previously completed and accepted independent 160-order BDF output, every 60s plus the same official execution endpoint; not every second and not claimed as newly run.'
    r['old_Q2_first3h_unrounded'],ro=cmp(v[:10801],ov,t[:10801])
    ow=old_workbook(repo/'q2_final_delivery/output/result2.xlsx');wcmp=rounded(v[1:10801])!=ow
    r['old_Q2_actual_workbook']={'rows_per_sheet':10800,'values_per_sheet':226800,'four_decimal_difference_count_TC':wcmp.sum(axis=(0,1)).tolist(),'max_abs_new_rounded_vs_old_cells_TC':abs(rounded(v[1:10801])-ow).max(axis=(0,1)).tolist()}
    q3_path=repo/'q3_refinement_delivery/output/solution.npz';q=np.load(q3_path);qt=q['time_s'];ix=np.searchsorted(t,qt);assert np.max(abs(t[ix]-qt))<1e-8
    r['old_Q3_timestamp_binary_max_difference_s']=float(np.max(abs(t[ix]-qt)))
    r['old_Q3_common60s_and_endpoint'],rq3=cmp(v[ix],q['profile_TC'],qt)
    r['paper_table_difference_count_TC']=[]
    for k,name in enumerate(['table_temperature.csv','table_moisture.csv']):
        tab=np.loadtxt(repo/'q2_final_delivery/output'/name,delimiter=',',skiprows=1)
        vals=rounded(v[(tab[:,0]*3600).astype(int),:,k][:,[0,5,10,15,20]])
        r['paper_table_difference_count_TC'].append(int(np.count_nonzero(vals!=tab[:,1:])))
    r['source_files']=[fp(p) for p in [new,reference,old_path,repo/'q2_final_delivery/output/result2.xlsx',q3_path]]
    r['interpretation']='All new T,C originate from one Appendix 3 cold-start solution. Differences from old Q2 reflect finite-volume Richardson versus primitive-flux collocation output, not parameter or physical-model changes. These are 1D numerical checks, not a 2D geometry error bound.'
    for name,a,b,tt,mask in [('reference_four_decimal_differences',v[ri],rv,rt,rr),('old_Q2_four_decimal_differences',v[:10801],ov,t[:10801],ro),('old_Q3_four_decimal_differences',v[ix],q['profile_TC'],qt,rq3)]:
        with (out/(name+'.csv')).open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['t_s','radius_cm','field','new_unrounded','comparison_unrounded','new_4dp','comparison_4dp'])
            for i,j,k in np.argwhere(mask):w.writerow([tt[i],j/10,['T_degC','C_dry_basis'][k],format(a[i,j,k],'.17g'),format(b[i,j,k],'.17g'),format(a[i,j,k],'.4f'),format(b[i,j,k],'.4f')])
    (out/'audit.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n');print(json.dumps(r,ensure_ascii=False,indent=2));return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1]/'inputs/historical_baseline');p.add_argument('--new',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.repo,a.new,a.reference,a.out)
