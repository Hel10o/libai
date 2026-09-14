#!/usr/bin/env python3
"""Extra receipt: compare the exact final decimal-formatting rule to independent references."""
from pathlib import Path
import argparse,json,numpy as np
from q1_solver import read_xlsx
from q1_validate import ramp_cylinder
from q1_excel import rounded,audit_excel

def main():
    b=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=b/'input'/'附件1.xlsx')
    p.add_argument('--out',type=Path,default=b/'output');a=p.parse_args()
    data,_=read_xlsx(a.input);z=np.load(a.out/'q1_unrounded.npz')
    audit_excel(a.out/'result1.xlsx',a.out/'q1_unrounded.npz',a.out)
    T=ramp_cylinder(data,z['time_s'],z['radius_cm']/100,800)
    C=np.load(a.out/'validation'/'independent_cheb_N240.npz')['C'];report={}
    for key,ref in [('temperature_degC',T),('moisture_dry_basis',C)]:
        x=np.array(rounded(z[key][1:]));y=np.array(rounded(ref[1:]))
        report[key]={'cells':int(x.size),'different_four_decimal_cells':int(np.count_nonzero(x!=y)),
                     'max_unrounded_difference':float(np.max(np.abs(z[key][1:]-ref[1:])))}
    assert all(x['different_four_decimal_cells']==0 for x in report.values())
    np.savez_compressed(a.out/'validation'/'production_thermal_Bessel800.npz',time_s=z['time_s'],radius_cm=z['radius_cm'],temperature_degC=T)
    (a.out/'rounding_reference_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
