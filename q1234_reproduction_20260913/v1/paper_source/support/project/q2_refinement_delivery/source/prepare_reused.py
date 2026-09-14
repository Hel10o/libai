"""Select unchanged inputs from the accepted handoff; never recompute a PDE."""
from pathlib import Path
import argparse,hashlib,json,shutil,zipfile
import numpy as np

def digest(p):
    b=Path(p).read_bytes()
    return {'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'git_blob_sha1':hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--handoff',type=Path,required=True);ap.add_argument('--arrays-zip',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    old=a.handoff/'q2_final_delivery';manifest=json.loads((old/'evidence/delivery_manifest.json').read_text('utf-8-sig'))
    entries={d['path']:d for d in manifest['files']}
    reduced={};sources=[]
    names=['baseline','hm_0.8','hm_1.2','hT_0.8','hT_1.2','pchip_all241','freeze_all_initial']
    with zipfile.ZipFile(a.arrays_zip) as z:
        for name in names:
            suffix=f'output/validation/sensitivity_{name}.npz'
            candidate=[x for x in z.namelist() if x.endswith(suffix)]
            if len(candidate)!=1:raise ValueError(f'Missing or ambiguous original array: {name}')
            b=z.read(candidate[0]);sha=hashlib.sha256(b).hexdigest()
            assert sha==entries[suffix]['sha256'],f'Original Pro array hash differs: {name}'
            import io
            with np.load(io.BytesIO(b),allow_pickle=False) as ar:
                assert np.array_equal(ar['time_s'],np.arange(10801))
                reduced['time_s']=ar['time_s'];reduced['radius_cm']=ar['radius_cm'][[0,10,20]]
                for f in ['temperature_degC','moisture_dry_basis']:
                    reduced[f'{name}__{f}']=ar[f][:,[0,10,20]]
                for f in ['average_temperature_degC','average_moisture_dry_basis']:
                    reduced[f'{name}__{f}']=ar[f]
                reduced[f'{name}__environment']=ar['environment']
            sources.append({'name':name,'original_zip_member':candidate[0],'sha256':sha,'bytes':len(b),'original_Pro_manifest_match':True,'kind':'original Pro, not local review recomputation'})
    np.savez_compressed(a.out/'sensitivity_tracks.npz',**reduced)
    report={'source_archive':a.arrays_zip.name,**digest(a.arrays_zip),'projection':'Exact selection: r=0,1,2 cm; stored rdr averages; environment; all seconds. No interpolation, rounding or PDE recomputation.','source_arrays':sources,'projected_file':digest(a.out/'sensitivity_tracks.npz')}
    (a.out/'sensitivity_provenance.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Selected and hash-verified',len(sources),'original Pro sensitivity arrays')
if __name__=='__main__':main()
