"""Q2 post-acceptance diagnostics, using existing unrounded trajectories only.

No PDE solver is called. Array identities and 1/2/5-second derivative stencils
are checked. Figures use no manually selected colors and each has its own axes.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,platform,sys,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

FIELDS=('temperature_degC','moisture_dry_basis')
POSITIONS=(0,10,20)
LABELS=('Axis, r = 0 cm','Interior, r = 1 cm','Surface, r = 2 cm')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def json_write(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

def environment_read(p):
    """Read raw worksheet values with stdlib ZIP/XML; no workbook modifications."""
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(p) as z:
        if z.testzip():raise ValueError('XLSX CRC failure')
        ss=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns):
                ss.append(''.join(i.text or '' for i in e.iter('{'+ns['s']+'}t')))
        wb=ET.fromstring(z.read('xl/workbook.xml'));sh=wb.find('s:sheets',ns)
        if len(sh)!=1:raise ValueError('Unexpected sheets')
        rid=sh[0].attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        target=next(x.attrib['Target'] for x in ET.fromstring(z.read('xl/_rels/workbook.xml.rels')) if x.attrib['Id']==rid)
        target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
        rows=[]
        for row in ET.fromstring(z.read(target)).findall('.//s:sheetData/s:row',ns):
            vals={}
            for c in row.findall('s:c',ns):
                if c.find('s:f',ns) is not None:raise ValueError('Formula in raw input')
                v=c.find('s:v',ns);typ=c.attrib.get('t')
                if typ=='s': val=ss[int(v.text)]
                elif typ=='inlineStr':val=''.join(e.text or '' for e in c.findall('.//s:t',ns))
                else:val=float(v.text) if v is not None else None
                vals[''.join(x for x in c.attrib['r'] if x.isalpha())]=val
            rows.append([vals.get(c) for c in 'ABC'])
    assert rows[0]==['时间','温度','水分浓度']
    a=np.asarray(rows[1:],float)
    assert a.shape==(241,3) and np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(0,14401,60))
    return a

def d_piecewise(y,h=1):
    """Fourth-order finite differences, not crossing the 60-s input knots.
    Central [-2h,-h,0,h,2h] away from knots, otherwise 5-point forward/backward.
    At a knot select the right segment; at the final endpoint select the left.
    Early 0..9 s values are not trusted as instantaneous-rate estimates.
    """
    n=y.shape[0];t=np.arange(n);last=n-1
    start=(t//60)*60;end=np.minimum(start+60,last)
    start[-1]=max(0,last-60)
    central=(t-2*h>=start)&(t+2*h<=end)
    forward=(~central)&(t+4*h<=end)
    backward=~central&~forward
    out=np.zeros_like(y)
    stencils=[(central,np.array([-2,-1,0,1,2]),np.array([1,-8,0,8,-1])/12),
              (forward,np.arange(5),np.array([-25,48,-36,16,-3])/12),
              (backward,-np.arange(5),-np.array([-25,48,-36,16,-3])/12)]
    for mask,steps,weights in stencils:
        idx=t[mask,None]+steps[None,:]*h
        # Subtract anchor before summation to reduce cancellation of constants.
        out[mask]=np.einsum('ijk,j->ik',y[idx]-y[t[mask],None,:],weights)/h
    return out

def run(out,raw,sens_path,env_path,validation_path,figures=True):
    out=Path(out);(out/'data').mkdir(parents=True,exist_ok=True);(out/'figures').mkdir(exist_ok=True)
    with np.load(raw,allow_pickle=False) as z:a={k:z[k] for k in z.files}
    t=a['time_s'];T=a[FIELDS[0]][:,POSITIONS];C=a[FIELDS[1]][:,POSITIONS];Tk=T+273.15
    assert np.array_equal(t,np.arange(10801)) and np.isfinite(T).all() and np.all(C>0)
    env=environment_read(env_path)
    np.testing.assert_allclose(a['environment'],np.column_stack([np.interp(t,env[:,0],env[:,j]) for j in (1,2)]),rtol=0,atol=0)
    np.savetxt(out/'data/environment_original.csv',env,delimiter=',',fmt='%.17g',header='time_s,temperature_degC,environmental_quantity_kgkg_basis_unspecified',comments='')
    D=.0024*np.exp(-.45/C)*np.exp(-3850/Tk);D0=float(.0024*np.exp(-.45/2.55)*np.exp(-3850/301.15))
    H=3850*(1/301.15-1/Tk);M=.45*(1/2.55-1/C);lnratio=np.log(D/D0)
    resid=float(np.max(abs(lnratio-H-M)))
    assert resid<2e-13
    derived={};ratecheck={}
    for h in [1,2,5]:
        td=d_piecewise(T,h);cd=d_piecewise(C,h)
        hd=3850/Tk**2*td;md=.45/C**2*cd
        ld=d_piecewise(lnratio,h)
        derived[h]=(td,cd,hd,md,hd+md,ld)
    trustworthy=(t>=10)&(t<=10790)
    # Additional conservative uncertainty band: all three step widths must agree.
    r1=derived[1][4];err=np.maximum(abs(r1-derived[2][4]),abs(r1-derived[5][4]))
    eps=np.maximum(5*err,1e-10)
    sign=np.where(r1>eps,1,np.where(r1<-eps,-1,0));sign[~trustworthy]=0
    ratecheck={'method':'piecewise fourth-order finite differences, 1/2/5 s stencils, no crossing 60 s knots; right derivative at knots; not a PDE RHS.',
      'early_rate_not_interpreted_s':[0,9],'t0_surface_rate_note':'Initial Robin corner produces a boundary layer; no finite instantaneous derivative is inferred from t=0.',
      'step_1_vs_2_max_abs_s-1_after10':float(np.max(abs(r1-derived[2][4])[trustworthy])),
      'step_1_vs_5_max_abs_s-1_after10':float(np.max(abs(r1-derived[5][4])[trustworthy])),
      'step_1_vs_5_max_abs_s-1_after60':float(np.max(abs(r1-derived[5][4])[t>=60])),
      'chain_rule_vs_direct_log_derivative_after10_s-1':float(np.max(abs(r1-derived[1][5])[trustworthy])),
      'sign_certification':'Pointwise |rate| > max(5*max(step differences),1e-10 s^-1). Diagnostic, not a formal error bound.'}
    columns=['time_s','radius_cm','temperature_degC','C_kgkg_dry','D_m2_s','H_log','M_log','D_over_D0','T_rate_K_s','C_rate_kgkg_s','H_rate_s-1','M_rate_s-1','lnD_rate_s-1','lnD_rate_step_error_s-1','net_sign_certified','rate_interpretable']
    rows=[]
    for j,p in enumerate(POSITIONS):
        rows.append(np.column_stack([t,np.full_like(t,p/10),T[:,j],C[:,j],D[:,j],H[:,j],M[:,j],D[:,j]/D0,
             *(derived[1][k][:,j] for k in range(5)),err[:,j],sign[:,j],trustworthy.astype(int)]))
    mat=np.concatenate(rows);np.savetxt(out/'data/mechanism_timeseries.csv',mat,delimiter=',',fmt='%.17g',header=','.join(columns),comments='')
    np.savez_compressed(out/'data/mechanism_unrounded.npz',time_s=t,radius_cm=np.array([0.,1.,2.]),T=T,C=C,D=D,H=H,M=M,D_over_D0=D/D0,
        rate_h1=r1,rate_h2=derived[2][4],rate_h5=derived[5][4],H_rate=derived[1][2],M_rate=derived[1][3],rate_step_difference=err,rate_interpretable=trustworthy)
    intervals=[];summary=[];table=[]
    for j,p in enumerate(POSITIONS):
        imax=int(np.argmax(D[:,j]));imin=int(np.argmin(D[:,j]));up=np.flatnonzero((lnratio[:-1,j]<0)&(lnratio[1:,j]>=0))
        back=float(t[up[0]]-lnratio[up[0],j]/(lnratio[up[0]+1,j]-lnratio[up[0],j])) if len(up) else None
        d={"radius_cm":p/10,'H_3h':float(H[-1,j]),'M_3h':float(M[-1,j]),'net_log_3h':float(lnratio[-1,j]),'D_3h':float(D[-1,j]),'D_ratio_3h':float(D[-1,j]/D0),
           'global_sample_min_s':imin,'D_ratio_min':float(D[imin,j]/D0),'global_sample_max_s':imax,'D_ratio_max':float(D[imax,j]/D0),'return_to_initial_D_s_linear_between_samples':back,
           'last_half_hour_delta_H':float(H[-1,j]-H[9000,j]),'last_half_hour_delta_M':float(M[-1,j]-M[9000,j]),'last_half_hour_D_relative_change':float(D[-1,j]/D[9000,j]-1),
           'rate_positive_time_fraction_2to3h':float(np.mean(sign[7200:10800,j]>0)),'rate_negative_time_fraction_2to3h':float(np.mean(sign[7200:10800,j]<0)),
           'temperature_cooling_seconds_2to3h':int(np.sum(derived[1][0][7200:10800,j]<-1e-9))}
        summary.append(d)
        for ti in range(0,10801,1800):table.append([ti/3600,p/10,T[ti,j],C[ti,j],H[ti,j],M[ti,j],D[ti,j],D[ti,j]/D0])
        for lo in range(0,10800,60):
            dh=H[lo+60,j]-H[lo,j];dm=M[lo+60,j]-M[lo,j]
            intervals.append([lo,lo+60,p/10,dh,dm,dh+dm,dh/60,dm/60,(dh+dm)/60])
    np.savetxt(out/'data/mechanism_table.csv',table,delimiter=',',fmt='%.17g',header='time_h,radius_cm,T_degC,C_kgkg,H,M,D_m2_s,D_over_D0',comments='')
    np.savetxt(out/'data/interval_60s_contributions.csv',intervals,delimiter=',',fmt='%.17g',header='start_s,end_s,radius_cm,delta_H,delta_M,delta_logD,mean_H_rate_s-1,mean_M_rate_s-1,mean_logD_rate_s-1',comments='')
    # Equilibrium mapping diagnostics keep the sampled state fixed; no extra PDE scenario.
    b=a['environment'][:,1];cs=C[:,2];q=8e-7*(cs-b);eb=-b/(cs-b)
    flux=np.column_stack([t,cs,b,cs-b,q,eb])
    np.savetxt(out/'data/boundary_identifiability_diagnostic.csv',flux,delimiter=',',fmt='%.17g',header='time_s,C_surface,b_equivalent,driving_gap,effective_flux_m_s,elasticity_to_b_at_fixed_state',comments='')
    # Reuse same-grid ±20% trajectories; these are finite secants, not fitted uncertainties.
    with np.load(sens_path,allow_pickle=False) as zz:s={k:zz[k] for k in zz.files}
    assert np.array_equal(s['time_s'],t)
    sr=[];elastic={};at_end=[]
    for par in ['hm','hT']:
        low,hi=par+'_0.8',par+'_1.2'
        Cb=np.column_stack([s['baseline__moisture_dry_basis'],s['baseline__average_moisture_dry_basis']])
        Tb=np.column_stack([s['baseline__temperature_degC'],s['baseline__average_temperature_degC']])
        Cl=np.column_stack([s[low+'__moisture_dry_basis'],s[low+'__average_moisture_dry_basis']])
        Ch=np.column_stack([s[hi+'__moisture_dry_basis'],s[hi+'__average_moisture_dry_basis']])
        Tl=np.column_stack([s[low+'__temperature_degC'],s[low+'__average_temperature_degC']])
        Th=np.column_stack([s[hi+'__temperature_degC'],s[hi+'__average_temperature_degC']])
        ec=(Ch-Cl)/(.4*Cb);et=np.divide(Th-Tl,.4*(Tb-28),out=np.full_like(Tb,np.nan),where=Tb-28>=.05)
        elastic[par]=ec
        for j,loc in enumerate(['axis','r1cm','surface','rdr_mean']):
            for ti in range(0,len(t)):
                sr.append([par,loc,int(ti),Cb[ti,j],Cl[ti,j],Ch[ti,j],ec[ti,j],Tb[ti,j],Tl[ti,j],Th[ti,j],et[ti,j] if np.isfinite(et[ti,j]) else ''])
            at_end.append({'parameter':par,'location':loc,'C_minus20':float(Cl[-1,j]),'C_baseline':float(Cb[-1,j]),'C_plus20':float(Ch[-1,j]),'C_secant_elasticity':float(ec[-1,j]),
              'T_minus20':float(Tl[-1,j]),'T_baseline':float(Tb[-1,j]),'T_plus20':float(Th[-1,j]),'T_rise_secant_elasticity':float(et[-1,j]),
              'maximum_abs_C_change':float(max(np.max(abs(Cl[:,j]-Cb[:,j])),np.max(abs(Ch[:,j]-Cb[:,j])))),
              'maximum_abs_T_change':float(max(np.max(abs(Tl[:,j]-Tb[:,j])),np.max(abs(Th[:,j]-Tb[:,j]))))})
    with (out/'data/parameter_sensitivity_timeseries.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['parameter','location','time_s','C_base','C_minus20','C_plus20','C_secant_elasticity','T_base','T_minus20','T_plus20','T_rise_secant_elasticity_blank_if_rise_lt_0.05K']);w.writerows(sr)
    with (out/'data/parameter_sensitivity_3h.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(at_end[0]));w.writeheader();w.writerows(at_end)
    # Empirical-density interpretation contradiction, not a new physical density estimate.
    rho=650+128*C;rho_d_implied=rho/(1+C);rho_d0=(650+128*2.55)/(1+2.55)
    volume=np.pi*.02**2*.25
    val=json.loads(Path(validation_path).read_text('utf-8-sig'))
    # Display only symbolic scaling for latent heat: energy = lambda * conditional lost mass.
    mloss_conditional=rho_d0*volume*(2.55-a['average_moisture_dry_basis'][-1])
    tail=env[env[:,0]>=10800]
    time_mean=np.trapezoid(tail[:,1:],tail[:,0],axis=0)/3600
    report={'scope':'New postprocessing only; no 3-hour PDE rerun and no corrected physical prediction.',
       'source_hashes':{'accepted_npz':sha(raw),'sensitivity_projection':sha(sens_path),'raw_environment':sha(env_path),'reused_validation':sha(validation_path)},
       'D_initial_m2_s':D0,'identity_max_abs_log_error':resid,'locations':summary,'derivative_check':ratecheck,'sensitivity_3h':at_end,
       'boundary_fixed_state':{'q0_effective_m_s':float(q[0]),'q3h_effective_m_s':float(q[-1]),'b_elasticity_initial':float(eb[0]),'b_elasticity_3h':float(eb[-1]),'no_boundary_mapping_PDE_scenarios':True},
       'density_diagnostic':{'interpretation':'Only if empirical rho is actual wet density, which is not the accepted static fixed-volume model.',
         'implied_rho_d_initial':rho_d0,'implied_rho_d_3h_0_1_2cm':rho_d_implied[-1].tolist(),
         'implied_ratio_to_initial_3h':(rho_d_implied[-1]/rho_d0).tolist(),
         'conditional_lost_water_kg_reference_dry_density':float(mloss_conditional),
         'reused_effective_heat_input_J':float(val['balance']['net_heat_input_J_m3']*volume),
         'latent_equal_effective_heat_input_at_lambda_J_kg':float(val['balance']['net_heat_input_J_m3']*volume/mloss_conditional),
         'no_lambda_assumed_and_no_latent_PDE_run':True},
       'time_scope':{'raw_data_s':[0,14400],'paper_s':[1800,3600,5400,7200,9000,10800],'excel_s':[1,10800],
        'tail_sample_mean_3to4h':tail[:,1:].mean(axis=0).tolist(),'tail_time_mean_linear_3to4h':time_mean.tolist(),'last_environment_value':env[-1,1:].tolist(),
        'no_days_long_solution_or_threshold_in_this_review':True},
       'response_3h':{'environment_T':float(a['environment'][-1,0]),'radial_T_span_K':float(np.ptp(T[-1])),
        'max_T_lag_to_current_air_K':float(np.max(abs(T[-1]-a['environment'][-1,0]))),
        'radial_C_span_kgkg':float(np.ptp(C[-1])),'remaining_effective_water_fraction':float(a['average_moisture_dry_basis'][-1]/2.55)},
       'runtime':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform()}}
    json_write(out/'data/mechanism_validation.json',report)
    if figures:draw(out,t,T,C,D,D0,H,M,derived,elastic,eb)
    print(json.dumps({'identity_error':resid,'locations':summary,'derivatives':ratecheck,'response_3h':report['response_3h']},ensure_ascii=False,indent=2))
    return report

def draw(out,t,T,C,D,D0,H,M,derived,elastic,eb):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.titlesize':13,'legend.fontsize':9,'svg.fonttype':'none'})
    def finish(fig,ax,name,title,ylabel):
        ax.set(xlabel='Time / h',ylabel=ylabel,title=title,xlim=(0,3));ax.grid(True,alpha=.25);ax.legend(loc='best');fig.tight_layout()
        fig.savefig(out/'figures'/f'{name}.png',dpi=300);fig.savefig(out/'figures'/f'{name}.svg');plt.close(fig)
    tx=t/3600
    fig,ax=plt.subplots(figsize=(8.2,4.7))
    for j in range(3):ax.plot(tx[::5],(D[:,j]/D0)[::5],label=LABELS[j])
    ax.axhline(1,ls=':',label='Initial diffusivity')
    finish(fig,ax,'01_D_competition','Local diffusivity: heating versus moisture loss','$D(t,r)/D_0$')
    for j in range(3):
        fig,ax=plt.subplots(figsize=(8.2,4.7));ax.plot(tx[::5],H[::5,j],label='Heating log contribution H')
        ax.plot(tx[::5],M[::5,j],label='Moisture log contribution M',ls='--');ax.plot(tx[::5],(H+M)[::5,j],label='Net ln(D / D0)',ls='-.')
        finish(fig,ax,f'0{j+2}_log_contributions','Cumulative contributions | '+LABELS[j],'Dimensionless log contribution')
    for j in range(3):
        fig,ax=plt.subplots(figsize=(8.2,4.7));mask=(t>=10)&(t<=10790)
        tt=tx[mask][::5]
        ax.plot(tt,derived[1][2][mask,j][::5]*3600,label='Thermal term (may become negative)')
        ax.plot(tt,derived[1][3][mask,j][::5]*3600,label='Moisture-loss term',ls='--')
        ax.plot(tt,derived[1][4][mask,j][::5]*3600,label='Net local rate',ls='-.')
        finish(fig,ax,f'0{j+5}_instantaneous_rates','Instantaneous-rate diagnostic | '+LABELS[j],'Rate of ln(D) / h$^{-1}$')
    for q,par in enumerate(['hm','hT']):
        fig,ax=plt.subplots(figsize=(8.2,4.7))
        for j,loc in enumerate(['Axis','r = 1 cm','Surface','r dr-weighted mean']):ax.plot(tx[::5],elastic[par][::5,j],label=loc)
        finish(fig,ax,f'0{q+8}_sensitivity_{par}',f'Moisture response to {par}: same-grid +/-20% scenarios','Finite-secant normalized sensitivity')
    fig,ax=plt.subplots(figsize=(8.2,4.7));ax.plot(tx[::5],eb[::5],label='b sensitivity at fixed simulated surface C');ax.axhline(0,ls=':')
    finish(fig,ax,'10_mapping_flux_diagnostic','Boundary-input sensitivity: algebraic diagnostic, not a new PDE','d ln(q) / d ln(b), holding surface C fixed')
    json_write(out/'figures/manifest.json',{'pairs':10,'format':['PNG 300 dpi','SVG'],'source':'q2_mechanism.py; accepted unrounded data only','files':[p.name for p in sorted((out/'figures').glob('*')) if p.suffix in ('.png','.svg')]})

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=root);p.add_argument('--raw',type=Path,default=root/'reused/q2_final_delivery/output/q2_unrounded.npz')
    p.add_argument('--sensitivity',type=Path,default=root/'reused/sensitivity_tracks.npz');p.add_argument('--environment',type=Path,default=root/'reused/q2_final_delivery/inputs/附件1.xlsx')
    p.add_argument('--validation',type=Path,default=root/'reused/q2_final_delivery/output/validation/validation.json');p.add_argument('--no-figures',action='store_true')
    a=p.parse_args();run(a.out,a.raw,a.sensitivity,a.environment,a.validation,not a.no_figures)
