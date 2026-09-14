#!/usr/bin/env python3
"""Generate figures, numerical tables and reviewable Chinese paper from run files."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from q1_solver import P,Environment,read_xlsx
from q1_excel import audit_excel
TIMES=[100,300,600,900,1200,1500,1800]

def md_table(headers,rows):
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---:']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)

def figures(out,data,raw,v):
    dest=out/'figures';dest.mkdir(exist_ok=True)
    listing=[]
    def save(name,title,xlabel,ylabel):
        ax=plt.gca();ax.set_title(title,fontsize=12,pad=12)
        ax.set_xlabel(xlabel,fontsize=11);ax.set_ylabel(ylabel,fontsize=11)
        ax.tick_params(labelsize=10);ax.grid(True,alpha=.22);ax.legend(fontsize=9)
        plt.tight_layout();plt.savefig(dest/(name+'.png'),dpi=300,bbox_inches='tight')
        plt.savefig(dest/(name+'.svg'),bbox_inches='tight');plt.close()
        listing.append(name)
    t=raw['time_s'];r=raw['radius_cm'];e=Environment(data);p=Environment(data,'pchip')
    sample=data[data[:,0]<=1800]
    for i,(name,title,label) in enumerate([
      ('01_environment_temperature','Drying-room temperature: measured boundary input','Ambient temperature (°C)'),
      ('02_environment_moisture','Drying-room moisture: equivalent boundary input','Equivalent concentration (kg/kg)')]):
        plt.figure(figsize=(7.2,4.5));plt.plot(t,e(t)[:,i],label='Piecewise linear (main)')
        plt.plot(t,p(t)[:,i],'--',label='PCHIP (sensitivity)')
        plt.plot(sample[:,0],sample[:,i+1],'o',ms=3,label='Attachment 1 measurements')
        save(name,title,'Time (s)',label)
    for key,field,label,name in [('temperature_degC','Temperature','Temperature (°C)','03_temperature_time'),
                                ('moisture_dry_basis','Dry-basis moisture','C (kg water / kg dry solid)','04_moisture_time')]:
        plt.figure(figsize=(7.2,4.5))
        for i in range(0,21,5):plt.plot(t,raw[key][:,i],label=f'r = {r[i]:g} cm')
        if key=='temperature_degC':plt.plot(t,raw['environment'][:,0],':',label='Drying room')
        save(name,field+' at the cylinder midplane','Time (s)',label)
    for f,label,name in [('T','Temperature (°C)','05_temperature_profiles'),('C','C (kg water / kg dry solid)','06_moisture_profiles')]:
        z=np.load(out/f'{f}_n10240.npz');rr=z['internal_radius_m']*100
        plt.figure(figsize=(7.2,4.5))
        for time,s in zip(z['snapshot_times_s'][1:],z['internal_snapshots'][1:]):
            plt.plot(rr,s,label=f't = {time:.0f} s')
        save(name,('Temperature' if f=='T' else 'Dry-basis moisture')+' profiles at the midplane','Radial distance from axis (cm)',label)
    g=v['grid_pairwise']+v['further_refinement']
    for f,unit,name in [('T','K','07_grid_temperature'),('C','kg/kg','08_grid_moisture')]:
        dx=np.array([100*P.R/x['n_coarse'] for x in g])
        full=np.array([x[f]['max_full'] for x in g]);tab=np.array([x[f]['max_table'] for x in g])
        plt.figure(figsize=(7.2,4.5));plt.loglog(dx,full,'o-',label='Full 1800 × 21 output: max |coarse - fine|')
        plt.loglog(dx,tab,'s--',label='Required 7 × 5 table: max |coarse - fine|')
        save(name,('Temperature' if f=='T' else 'Moisture')+' grid refinement','Coarse radial spacing (cm)',f'Successive-grid difference ({unit})')
    (dest/'figure_manifest.json').write_text(json.dumps({'files':listing,'formats':['300 dpi PNG','SVG'],'source':'actual unrounded production and validation outputs'},indent=2),encoding='utf-8')
    return listing

def main():
    b=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,default=b/'input'/'附件1.xlsx')
    ap.add_argument('--out',type=Path,default=b/'output');args=ap.parse_args();out=args.out
    z=dict(np.load(out/'q1_unrounded.npz'));v=json.loads((out/'validation'/'validation.json').read_text())
    run=json.loads((out/'run.json').read_text());data,audit=read_xlsx(args.input)
    ex=audit_excel(out/'result1.xlsx',out/'q1_unrounded.npz',out)
    figures(out,data,z,v)
    tokens={}
    for key,title in [('temperature_degC','TEMP_TABLE'),('moisture_dry_basis','MOIST_TABLE')]:
        tokens[title]=md_table(['时间/s','0 cm','0.5 cm','1 cm','1.5 cm','2 cm'],[[str(t)]+[f'{x:.4f}' for x in z[key][t,::5]] for t in TIMES])
    g=v['grid_pairwise']+v['further_refinement']
    tokens['GRID_TABLE']=md_table(['区间数 N→2N','粗网格间距/cm','全场温差/K','全场含水率差/(kg/kg)','两张论文表最大含水率差'],
        [[f"{x['n_coarse']}→{x['n_fine']}",f"{100*P.R/x['n_coarse']:.7g}",f"{x['T']['max_full']:.3e}",f"{x['C']['max_full']:.3e}",f"{x['C']['max_table']:.3e}"] for x in g])
    indT=v['production_actual_thermal_analytic'];indC=v['independent_nonlinear_collocation']['fv10240_vs_N240']
    tests=v['tests'];lt=tests['latent_conditional'];sc=v['scales'];two=v['2d_heat_refinement']
    vals={
      'AVG_T':z['average_temperature_degC'][-1],'AVG_C':z['average_C'][-1],
      'CENTER_T':z['temperature_degC'][-1,0],'SURFACE_T':z['temperature_degC'][-1,-1],
      'CENTER_C':z['moisture_dry_basis'][-1,0],'SURFACE_C':z['moisture_dry_basis'][-1,-1],
      'LOSS_PERCENT':100*(P.C0-z['average_C'][-1])/P.C0,'TEMP_DIFF':z['temperature_degC'][-1,-1]-z['temperature_degC'][-1,0],
      'MEAN_2D_T':two['average_2d_1800'],'END_T':two['end_center_1800'],
      'LATENT_RATIO':lt['ratio_latent_to_sensible'],'LATENT_LOSS_KG':lt['implied_water_loss_kg'],
      'LATENT_J':lt['implied_latent_J'],'SENSIBLE_J':lt['baseline_sensible_J'],
      'LATENT_INITIAL_FLUX':lt['initial_latent_W_m2'],'LATENT_T_CENTER':lt['T1800_degC'][0],
      'LATENT_T_SURFACE':lt['T1800_degC'][-1],
      'HM_LOW':tests['hm_factor_0.8']['C1800'][-1],'HM_HIGH':tests['hm_factor_1.2']['C1800'][-1],
    }
    tokens.update({k:f'{float(x):.4f}' for k,x in vals.items()})
    for key,val in {
      'TIME_T':v['production_time_precision']['T']['max_full'],'TIME_C':v['production_time_precision']['C']['max_full'],
      'INDEP_T':indT['max_full'],'INDEP_C':indC['max_full'],'INDEP_C_TABLE':indC['max_table'],
      'CONSERV_T':run['T']['conservation_augmented_max_abs'],'CONSERV_C':run['C']['conservation_augmented_max_abs'],
      'QUAD_T':run['T']['independent_flux_quadrature_max_abs'],'QUAD_C':run['C']['independent_flux_quadrature_max_abs'],
      'NOFLUX_T':tests['no_flux_T']['integral_average_drift'],'NOFLUX_C':tests['no_flux_C']['integral_average_drift'],
      'PCHIP_T':tests['pchip_T']['max_full'],'PCHIP_C':tests['pchip_C']['max_full'],
      'TWO_D_DIFF':two['400_midplane_vs_production_K'],
      'TWO_D_TRUNC':two['200_vs_400_all_max_K'],
      'CONST_T':tests['constant_analytic_T']['fv_vs_series']['max_full'],
      'CONST_C':tests['constant_analytic_C']['fv_vs_series']['max_full'],
    }.items():tokens[key]=f'{val:.4e}'
    tokens['TIMESTEPS_T']=str(run['T']['accepted_steps']);tokens['TIMESTEPS_C']=str(run['C']['accepted_steps'])
    tokens['VERSIONS']=', '.join(k+' '+s for k,s in run['versions'].items())
    template=(Path(__file__).parent/'paper_template.md').read_text(encoding='utf-8')
    for k,val in tokens.items():template=template.replace('@@'+k+'@@',val)
    if '@@' in template:raise ValueError('Unresolved report token')
    paper=out/'第一问论文正文.md';paper.write_text(template,encoding='utf-8')
    # Audit the two Markdown tables against all corresponding workbook cells.
    for key,name in [('temperature_degC','TEMP_TABLE'),('moisture_dry_basis','MOIST_TABLE')]:
        assert tokens[name] in template
    summary={'main_model':'Midplane 1D effective Robin, independent sensible-heat and nonlinear dry-basis diffusion',
      'production_grid_intervals':10240,'values_1800':{k:float(x) for k,x in vals.items()},
      'independent_comparison':{'T':indT,'C':indC},'excel':ex['checks'],
      'paper_tables_generated_from_unrounded_output':True,'paper_tables_same_as_Excel':True,
      'limitations':['Original problem PDF formula image not read; repository extracted text and explicit supplied formula used',
        'Foods11 LFS full PDF not obtained','No internal-response experiment','No 2D nonlinear moisture computation',
        'Latent omission not small under literal evaporating-mass interpretation; baseline is conditional']}
    (out/'delivery_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Generated:',paper,'and 8 PNG/SVG figure pairs',flush=True)
if __name__=='__main__':main()
