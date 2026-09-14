"""Build six validated tables and five figures from frozen evidence only.

No PDE solver is imported or executed. Writes are restricted to the asset
directories and three named evidence files in this draft.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, FixedFormatter, MultipleLocator
from openpyxl import load_workbook

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TABLES=HERE/"tables"
FIGURES=HERE/"figures"
EVIDENCE=HERE/"evidence"
COLORS=["#0072B2","#D55E00","#009E73","#CC79A7","#666666"]
STYLES=["-","--","-.",":",(0,(5,1,1,1))]
SOURCES={}
OUTPUTS=[]
MESSAGES=[]
EXPECTED={}
FIGURE_RECORDS=[]


def log(message):
    print(message,flush=True)
    MESSAGES.append(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path):
    return path.relative_to(ROOT).as_posix()


def record_source(path, purpose):
    key=rel(path)
    if key not in SOURCES:SOURCES[key]={"sha256":sha(path),"purposes":[]}
    if purpose not in SOURCES[key]["purposes"]:SOURCES[key]["purposes"].append(purpose)
    return path


def load_npz(path,purpose):
    record_source(path,purpose)
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}


def csv_rows(path,purpose):
    record_source(path,purpose)
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.reader(f))


def match_time(times,wanted):
    i=int(np.argmin(abs(times-wanted)))
    assert abs(float(times[i])-float(wanted))<1e-6,(wanted,times[i])
    return i


def fmt(value):
    return "--" if value is None or not np.isfinite(value) else f"{float(value):.4f}"


def make_table(name,times,values,time_unit,source_files,official_rows):
    body=[]
    for t,values_row in zip(times,values):
        time_text=str(int(t)) if time_unit=="s" else fmt(t)
        body.append([time_text,*[fmt(v) for v in values_row]])
    # Official CSVs are independent frozen table artifacts. Compare display
    # values rather than relying on agreement of underlying float encodings.
    assert len(body)==len(official_rows)
    for actual,official in zip(body,official_rows):
        expected_time=str(int(float(official[0]))) if time_unit=="s" else fmt(float(official[0]))
        expected=[expected_time,*[fmt(float(v)) if v.strip() else "--" for v in official[1:1+len(actual)-1]]]
        assert actual==expected,(name,actual,expected)
    q4=name=="q4_moisture"
    columns=7 if q4 else 6
    lines=[r"\begin{tabular}{"+"r"*columns+"}",r"\toprule"]
    if q4:
        lines += [r"时间 / h & \multicolumn{5}{c}{径向距离 / cm} & 当前表面 \\",
                  r"\cmidrule(lr){2-6}",r" & 0 & 0.5 & 1 & 1.5 & 2 & \\"]
    else:
        lines += [f"时间 / {time_unit}"+r" & \multicolumn{5}{c}{径向距离 / cm} \\",
                  r"\cmidrule(lr){2-6}",r" & 0 & 0.5 & 1 & 1.5 & 2 \\"]
    lines += [r"\midrule",*[" & ".join(row)+r" \\" for row in body],r"\bottomrule",r"\end{tabular}"]
    path=TABLES/(name+".tex")
    path.write_text("\n".join(lines)+"\n",encoding="utf-8")
    OUTPUTS.append(path)
    EXPECTED[name]={"latex_file":rel(path),"column_count":columns,
        "headers":["时间/"+time_unit,"0 cm","0.5 cm","1 cm","1.5 cm","2 cm"]+(["当前表面"] if q4 else []),
        "time_unit":time_unit,"rows":body,
        "flat_numeric_cells":[value for row in body for value in row if value!="--"],
        "body_numeric_cell_count":sum(value!="--" for row in body for value in row),
        "moisture_or_temperature_numeric_cell_count":sum(value!="--" for row in body for value in row[1:]),
        "domain_empty_cell_count":sum(value=="--" for row in body for value in row),
        "sources":[rel(p) for p in source_files],"all_frozen_csv_four_decimal_cells_match":True}
    log(f"TABLE PASS {name}: {len(body)} rows, {columns} columns, {EXPECTED[name]['domain_empty_cell_count']} domain-empty cells")


def style():
    available={item.name for item in font_manager.fontManager.ttflist}
    selected=next((s for s in ("Microsoft YaHei","SimSun") if s in available),None)
    if selected is None:raise RuntimeError("Chinese plotting font not found")
    plt.rcParams.update({"font.family":selected,"font.size":8.5,"axes.labelsize":8.5,
        "axes.titlesize":9,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":7.5,
        "axes.linewidth":.65,"lines.linewidth":1.3,"axes.unicode_minus":False,
        "pdf.fonttype":42,"ps.fonttype":42,"figure.facecolor":"white","savefig.dpi":300})
    return selected


def panels(width=16.,height=7.4):
    fig,axes=plt.subplots(1,2,figsize=(width/2.54,height/2.54))
    fig.subplots_adjust(left=.085,right=.975,bottom=.19,top=.87,wspace=.32)
    for ax in axes:
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",color="#D7DCE0",lw=.5,alpha=.8)
        ax.set_axisbelow(True)
        ax.tick_params(direction="out",pad=2)
    return fig,axes


def save_figure(fig,name,record):
    for suffix in ("pdf","png"):
        path=FIGURES/(name+"."+suffix)
        if suffix=="pdf":fig.savefig(path,metadata={"CreationDate":None,"ModDate":None})
        else:fig.savefig(path,dpi=300)
        OUTPUTS.append(path)
    plt.close(fig)
    FIGURE_RECORDS.append({"name":name,**record})
    log("FIGURE BUILT "+name)


def read_radius():
    path=record_source(ROOT/"q4_complete_delivery/v1/inputs/attachment2.xlsx","q4_shrinkage")
    wb=load_workbook(path,read_only=True,data_only=True)
    data=np.asarray([row[:2] for row in list(wb.active.values)[1:] if row[0] is not None],float)
    wb.close()
    assert data.shape==(145,2) and np.all(np.diff(data[:,0])>0)
    return path,data


def main():
    for directory in (TABLES,FIGURES,EVIDENCE):directory.mkdir(parents=True,exist_ok=True)
    font=style()
    log("POSTPROCESS ONLY: no PDE solver imported or run")
    log(f"RUNTIME Python={platform.python_version()} numpy={np.__version__} scipy={scipy.__version__} matplotlib={matplotlib.__version__}; font={font}; executable={sys.executable}")
    q1path=ROOT/"q1_complete_delivery/q1_delivery/output/q1_unrounded.npz"
    q23path=ROOT/"q4_complete_delivery/v1/q123_closeout/output/q23_unified.npz"
    q3path=ROOT/"q3_refinement_delivery/output/solution.npz"
    q4path=ROOT/"q4_complete_delivery/v1/output/main.npz"
    q1=load_npz(q1path,"Q1 actual table states")
    q23=load_npz(q23path,"Current Q2 table states and Q3 full-process plot")
    q3=load_npz(q3path,"Current formal Q3 table 5 and endpoint crosscheck")
    q4=load_npz(q4path,"Current Q4 table 6, full-process and radial plots")
    positions=[0,5,10,15,20]
    assert np.array_equal(q1["radius_cm"][positions],np.array([0,.5,1,1.5,2]))
    assert np.array_equal(q23["radius_cm"][positions],np.array([0,.5,1,1.5,2]))
    for suffix,key in (("temperature","temperature_degC"),("moisture","moisture_dry_basis")):
        official_path=ROOT/f"q1_complete_delivery/q1_delivery/output/table_{suffix}.csv"
        official=csv_rows(official_path,"Q1 frozen official four-decimal crosscheck")[1:]
        times=np.array([100,300,600,900,1200,1500,1800])
        values=q1[key][[match_time(q1["time_s"],t) for t in times]][:,positions]
        make_table("q1_"+suffix,times,values,"s",[q1path,official_path],official)
    for suffix,channel in (("temperature",0),("moisture",1)):
        official_path=ROOT/f"q2_final_delivery/output/table_{suffix}.csv"
        official=csv_rows(official_path,"Q2 frozen official four-decimal crosscheck")[1:]
        times=np.arange(.5,3.01,.5)
        values=q23["profile_TC"][[match_time(q23["time_s"],t*3600) for t in times]][:,positions,channel]
        make_table("q2_"+suffix,times,values,"h",[q23path,official_path],official)
    q3csv=ROOT/"q3_refinement_delivery/output/table5.csv"
    official=csv_rows(q3csv,"Current formal Q3 table 5 display crosscheck")[1:]
    times=np.array([float(row[0]) for row in official])
    values=q3["profile_TC"][[match_time(q3["time_s"],t*3600) for t in times]][:,positions,1]
    make_table("q3_moisture",times,values,"h",[q3path,q3csv],official)
    # The new unified solution gives exactly the same four-decimal Q3 rows.
    common=q23["profile_TC"][[match_time(q23["time_s"],t*3600) for t in times]][:,positions,1]
    assert [[fmt(v) for v in row] for row in common]==[[fmt(v) for v in row] for row in values]
    log("Q3 current formal table and new Q23 unified trajectory: all four-decimal cells agree")
    q4csv=ROOT/"q4_complete_delivery/v1/output/table6_unrounded.csv"
    official=csv_rows(q4csv,"Current formal Q4 table 6 display crosscheck")[1:]
    times=np.array([float(row[0]) for row in official])
    indices=[match_time(q4["time_s"],t*3600) for t in times]
    values=np.c_[q4["sample_TC"][indices][:,positions,1],q4["surface_TC"][indices,1]]
    assert np.array_equal(np.isnan(values[:,:5]),q4["fixed_radius_m"][positions][None,:]>q4["radius_m"][indices,None]+1e-12)
    make_table("q4_moisture",times,values,"h",[q4path,q4csv],official)
    EXPECTED["q4_moisture"]["row_actual_radius_cm"]=[fmt(v*100) for v in q4["radius_m"][indices]]

    for name in ("q1_profiles","q2_profiles"):
        copied=[]
        for suffix in ("pdf","png"):
            original=record_source(ROOT/f"paper/q1_q2_stage/figures/{name}.{suffix}","Accepted profile figure copied byte-for-byte")
            target=FIGURES/original.name
            shutil.copyfile(original,target)
            assert sha(original)==sha(target)
            OUTPUTS.append(target)
            copied.append({"from":rel(original),"to":rel(target),"sha256":sha(original)})
        FIGURE_RECORDS.append({"name":name,"operation":"byte-for-byte copy of accepted figure","copied":copied})
        log("FIGURE COPIED "+name)

    fig,axes=panels(height=7.5)
    qt=np.unique(np.r_[np.arange(0,206907,60),206906.76])
    qi=np.array([match_time(q23["time_s"],t) for t in qt])
    series=[("Q3：固定半径",qt/3600,q23["profile_TC"][qi,0,1],q23["profile_TC"][qi,-1,1],57.4741),
            ("Q4：规定收缩",q4["time_s"]/3600,q4["full_TC"][:,0,1],q4["surface_TC"][:,1],51.0921)]
    for ax,(title,t,center,surface,execution) in zip(axes,series):
        ax.plot(t,center,color=COLORS[0],label="轴心")
        ax.plot(t,surface,color=COLORS[1],ls="--",label="当前表面")
        ax.axhline(.15,color="#666666",lw=.9,ls=":",label="达标阈值 0.15")
        ax.plot(execution,center[-1],"o",color=COLORS[0],ms=3.5)
        ax.annotate(f"{execution:.4f} h",xy=(execution,center[-1]),xytext=(execution-14,.245),
            fontsize=7.5,arrowprops={"arrowstyle":"-","lw":.6,"color":"#555555"})
        ax.set(xlabel="时间 / h",ylabel="干基含水率 / (kg/kg)",title=title,xlim=(0,60),ylim=(.045,2.9),yscale="log")
        ticks=[.05,.1,.15,.5,1.,2.5]
        ax.yaxis.set_major_locator(FixedLocator(ticks));ax.yaxis.set_major_formatter(FixedFormatter(["0.05","0.1","0.15","0.5","1","2.5"]))
        ax.minorticks_off();ax.xaxis.set_major_locator(MultipleLocator(12))
        ax.legend(loc="upper right",frameon=False,handlelength=2.2)
    save_figure(fig,"q3_q4_drying",{"sources":[rel(q23path),rel(q4path)],
        "panels":"Q3 left; Q4 right; shared x 0-60 h; logarithmic C axis; axis/surface/0.15 threshold; curves stop at official execution",
        "time_sampling":"Q3 every 60 s plus execution; Q4 all stored 60 s outputs plus execution"})

    radius_path,observed=read_radius()
    fig,axes=panels(height=7.5)
    axes[0].plot(observed[:,0]/3600,observed[:,1],color=COLORS[0],lw=1.2,label="附件2观测")
    axes[0].scatter(observed[::4,0]/3600,observed[::4,1],s=5,color=COLORS[0],zorder=3)
    axes[0].axvline(51.0921,color=COLORS[1],ls="--",lw=1,label="执行 51.0921 h")
    axes[0].set(xlabel="时间 / h",ylabel="外半径 / cm",title="(a) 观测半径历程",xlim=(0,72),ylim=(1.15,2.04))
    axes[0].xaxis.set_major_locator(MultipleLocator(12));axes[0].legend(frameon=False,loc="upper right")
    profile_hours=[6,12,24,51.0921]
    for h,color,ls in zip(profile_hours,COLORS,STYLES):
        i=match_time(q4["time_s"],h*3600)
        r_cm=q4["radius_m"][i]*100*np.sqrt(q4["x"])
        c=q4["full_TC"][i,:,1]
        axes[1].plot(r_cm,c,color=color,ls=ls,label=(f"{h:g} h" if h!=51.0921 else "结束 51.0921 h"))
        axes[1].plot(r_cm[-1],c[-1],"o",ms=2.5,color=color)
        assert np.all(r_cm<=q4["radius_m"][i]*100+1e-12)
    axes[1].axhline(.15,color="#aaaaaa",lw=.7,ls=":")
    axes[1].set(xlabel="当前径向距离 / cm",ylabel="干基含水率 / (kg/kg)",title="(b) 当前材料域内剖面",xlim=(0,1.45),ylim=(0,1.84))
    axes[1].xaxis.set_major_locator(MultipleLocator(.5));axes[1].legend(frameon=False,loc="upper right",fontsize=7.)
    save_figure(fig,"q4_shrinkage",{"sources":[rel(radius_path),rel(q4path)],"profile_hours":profile_hours,
        "panels":"left observed R over 0-72 h; right actual-r profiles, points terminate at current surface; no extrapolation outside material"})

    fig,axes=panels(height=7.5)
    scenarios=["temperature_minus_1K","temperature_plus_1K","boundary_minus_0p005","boundary_plus_0p005"]
    labels=["温度 −1 K","温度 +1 K","边界 b −0.005","边界 b +0.005"]
    plotdata=[]
    for ax,case,title in zip(axes,("q3","q4"),("Q3：固定半径","Q4：规定收缩")):
        base_path=record_source(ROOT/f"q1234_overall_review_delivery/v1/evidence/environment/{case}_baseline_n40.json","Same-n40 environment baseline")
        baseline=json.loads(base_path.read_text(encoding="utf-8"))
        changes=[]
        for scenario in scenarios:
            path=record_source(ROOT/f"q1234_overall_review_delivery/v1/evidence/environment/{case}_{scenario}_n40.json","Executed environment root difference")
            data=json.loads(path.read_text(encoding="utf-8"))
            assert data["case"]==case and data["n"]==40
            change=(data["critical_s"]-baseline["critical_s"])/3600
            changes.append(change);plotdata.append({"case":case,"scenario":scenario,"delta_root_h":change,"source":rel(path)})
        y=np.arange(4)
        ax.barh(y,changes,color=[COLORS[1] if d>0 else COLORS[0] for d in changes],height=.56)
        ax.axvline(0,color="#666666",lw=.8)
        for yy,value in zip(y,changes):
            if abs(value)>.3:
                # Large-bar labels sit inside bars, avoiding overlap with the
                # Chinese y tick labels outside the plotting rectangle.
                ax.text(value-(.07 if value>=0 else -.07),yy,f"{value:+.3f}",va="center",ha="right" if value>=0 else "left",fontsize=7.3,color="white")
            else:
                ax.text(value+(.06 if value>=0 else -.06),yy,f"{value:+.3f}",va="center",ha="left" if value>=0 else "right",fontsize=7.3)
        ax.set_yticks(y,labels,fontsize=7.3);ax.invert_yaxis()
        ax.set(xlabel="临界根变化 / h",title=title,xlim=(-2.45,2.45))
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.grid(axis="y",visible=False);ax.grid(axis="x",color="#D7DCE0",lw=.5)
    fig.subplots_adjust(left=.16,right=.975,wspace=.60)
    save_figure(fig,"environment_sensitivity",{"panels":"Both are analyst-selected changes after 4 h, not error bars or confidence intervals; Q4 R(t) unchanged", "data":plotdata})

    expected_path=EVIDENCE/"table_expected.json"
    expected_document={"schema_version":1,"scope":"Body cells, including times; all nonblank values are exact desired PDF strings.",
        "tables":EXPECTED,"total_numeric_body_cells":sum(v["body_numeric_cell_count"] for v in EXPECTED.values()),
        "total_temperature_moisture_cells":sum(v["moisture_or_temperature_numeric_cell_count"] for v in EXPECTED.values()),
        "total_domain_empty_cells":sum(v["domain_empty_cell_count"] for v in EXPECTED.values())}
    expected_path.write_text(json.dumps(expected_document,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    OUTPUTS.append(expected_path)
    for source,description in SOURCES.items():assert sha(ROOT/source)==description["sha256"],source
    assert len(EXPECTED)==6 and len(FIGURE_RECORDS)==5
    log(f"FINAL PASS: 6 tables; {expected_document['total_temperature_moisture_cells']} temperature/moisture cells and {expected_document['total_domain_empty_cells']} domain-empty cells; 5 figure PDFs + PNG previews; all frozen source hashes unchanged")
    logfile=EVIDENCE/"assets_build.log"
    logfile.write_text("\n".join(MESSAGES)+"\n",encoding="utf-8")
    OUTPUTS.append(logfile)
    manifest={"git_commit_read":"anonymous-code-snapshot",
        "script":rel(Path(__file__)),"script_sha256":sha(Path(__file__)),
        "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__,"executable":sys.executable,"font":font},
        "scope":"Frozen evidence postprocessing only; no PDE runs; generated files replace only named draft asset outputs.",
        "sources":SOURCES,"figures":FIGURE_RECORDS,
        "outputs":{rel(path):{"sha256":sha(path),"bytes":path.stat().st_size} for path in OUTPUTS}}
    (EVIDENCE/"assets_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":main()
