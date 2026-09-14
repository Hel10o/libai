"""Build geometry discussion from completed comparison JSON and endpoint metadata."""
from pathlib import Path
import json
import numpy as np
from geometry_solver import DEST
from geometry_analysis import summarize

rows,data=summarize();v1=DEST.parents[1]
q3=json.loads((v1/'q123_closeout/validation/q2_closeout_summary.json').read_text())
q4=json.loads((v1/'output/end_event.json').read_text())
main4=json.loads((v1/'output/main.json').read_text())
def val(key,group='full_common_range'):
    return data[key]['intervals'][group]
def table(selected):
    out=['|模型/通量|径向×轴向区间|配对一维临界/s|二维临界/s|二维−配对一维/s|',
         '|---|---:|---:|---:|---:|']
    for key in selected:
        if key not in data:continue
        s=data[key]['summary'];tag=s['case']+'/'+s['flux']+('（紧配置）' if s['rtol']!=2e-9 else '')
        out.append(f"|{tag}|{s['nr']}×{s['nz']}|{s['event_1d_s']:.9f}|{s['event_2d_s']:.9f}|{s['end_effect_s']:.9f}|")
    return '\n'.join(out)

q3keys=['q23_40x64','q23_80x128_integral','q23_80x256_integral']
q4keys=['q4_40x64','q4_40x64_integral','q4_40x128_integral','q4_80x128_integral','q4_40x64_integral_rtol2e-11']
q3d=[data[k]['summary']['end_effect_s'] for k in q3keys[1:] if k in data]
q4d=[data[k]['summary']['end_effect_s'] for k in q4keys[1:] if k in data]
q3est=[q3['critical_s']+min(q3d),q3['critical_s']+max(q3d)]
q4est=[q4['critical_s']+min(q4d),q4['critical_s']+max(q4d)]
# Explicit estimate: translate the one-dimensional local event slope by the
# matched geometry shift; these are not evaluated coarse 2D endpoint fields.
slope=main4['stats']['critical_slope'] if 'stats' in main4 and 'critical_slope' in main4['stats'] else None
if slope is None:
    def find(x):
        if isinstance(x,dict):
            if 'critical_slope' in x:return x['critical_slope']
            for v in x.values():
                a=find(v)
                if a is not None:return a
        return None
    slope=find(main4)
c4est=sorted(.15+slope*(q4['execution_s']-t) for t in q4est)
est={'method':'matched 2D-minus-1D event correction applied to converged 1D; Cmax uses local event linearization. Estimates, not full-field evaluations at official execution time.',
     'q3_converged_1d_root_s':q3['critical_s'],'q3_corrected_2d_root_range_s':q3est,
     'q4_converged_1d_root_s':q4['critical_s'],'q4_matched_geometry_shift_range_s':[min(q4d),max(q4d)],
     'q4_corrected_2d_root_range_s':q4est,'q4_execution_s':q4['execution_s'],
     'q4_linearized_Cmax_estimate_at_execution':c4est,'q4_1d_critical_slope':slope,
     'all_completed_q3_q4_event_wettest_locations_r_z_m':{k:data[k]['end_location_r_z_m'] for k in q3keys+q4keys if k in data}}
(DEST/'corrected_global_event_estimates.json').write_text(json.dumps(est,indent=2,ensure_ascii=False)+'\n')

a=val('q1_80x128');b=val('q23_80x128_integral','q2_table_range')
finite_q3='q23_80x256_integral' if 'q23_80x256_integral' in data else 'q23_80x128_integral'
new_q3=data[finite_q3]['summary']
paragraph=f'''# 同一有效模型的一维/二维适用范围与第三问终点

本轮固定读取提交为 `anonymous-source`。比较均采用用户选择的题给物性、有效热容量与等效空气含水率边界，不含显式潜热；没有把 F 或 R_L 的旧二维差混入当前模型的误差依据。

## 方法与可比较性

二维区域为圆柱的一半轴向剖面，`0≤r≤R`、`0≤z≤0.125 m`；`z=0`为中截面。轴心及中截面施加对称条件，侧面和端面均使用 `hT=25`、`hm=8e-7`，角部按真实相交面面积求和，没有新增端面倍率。

每个二维算例都从统一初态 `T=28 ℃、C=2.55` 联立求温湿场，并配对完全相同径向节点、径向通量和积分配置的一维算例。由“二维临界−同网格一维临界”提取端部效应，避免把粗径向离散偏差误称为几何影响。

节点对偶有限体积采用二次加密网格，向侧面和端面加密；热通量使用局部导热系数的调和面值，主要水通量使用 `D=A(T) exp(-B/C)` 的含水率积分。另保留调和扩散系数的独立粗离散作为交叉证据。

4h 前在附件1每个60s节点分段，之后为末小时61点算术均值平台。全程温度都参与联立更新，没有后期等温投影。以下最大差是每60s和题目指定时刻、21个径向输出点的采样差，不能冒充连续时空严格上界。

## 第一、二问：中截面与端部应分别说明

|范围与本轮网格|中截面最大温差/K|中截面最大含水率差/(kg/kg)|端面最大温差/K|端面最大含水率差/(kg/kg)|
|---|---:|---:|---:|---:|
|Q1，0–1800s，80×128|{a['max_mid_difference_T_C'][0]:.8g}|{a['max_mid_difference_T_C'][1]:.8g}|{a['max_end_difference_T_C'][0]:.8g}|{a['max_end_difference_T_C'][1]:.8g}|
|Q2，0–10800s，80×128|{b['max_mid_difference_T_C'][0]:.8g}|{b['max_mid_difference_T_C'][1]:.8g}|{b['max_end_difference_T_C'][0]:.8g}|{b['max_end_difference_T_C'][1]:.8g}|

Q1 的40×64与80×128配对中截面温差分别约 `8.46e-8` 与 `5.65e-8 K`，含水率差均约 `4e-9`，已接近积分误差水平。其一维解足以代表题表指定的中截面径向输出，结论不推广到端面或整体平均值。

Q2 中截面温差约 `0.00169 K`，因此不能声称“一维与二维所有温度均四位小数一致”。含水率差约 `2.25e-5 kg/kg`；接近舍入边界的单格仍可能变化。工程建模容差、四位显示精度、离散误差与几何差属于不同尺度。

## 第三问：有限圆柱的全域临界

{table(q3keys)}

原无潜热一维模型的精细临界为 `{q3['critical_s']:.9f} s`；本轮第二问完整逐秒轨迹继续导出到 `{q3['time_end_s']:.2f} s`，对应 `57.4741 h`。这一端点仍属该有效模型的保守执行口径，不把它重新命名为精确二维最短时长。

最细已完成二维 `{new_q3['nr']}×{new_q3['nz']}` 的原始临界是 `{new_q3['event_2d_s']:.9f} s`，最湿节点确为 `(r,z)=(0,0)`。同一原始场在等号根前后1s分别高于与低于阈值，完整场保存在相应 `event_states` 数组。

去除匹配的一维径向偏差后，当前已完成积分网格给出的有限圆柱临界估计范围为 `{q3est[0]:.6f}–{q3est[1]:.6f} s`。这些数是相同模型的配对校正估计，不是另一次精细二维实跑，也不是严格误差界；它们支持保留原一维执行端点的保守方向。

历史同模型40×64、80×64、40×128结果已检查：端部效应约 `-7.66`、`-7.66`、`-7.51 s`，其中后两组含已量化的微小后期等温投影。本轮全程耦合80×128补齐交叉节点并保存了完整二维事件场；历史投影结果只用于核对趋势。

本轮拟做的160×128可选交叉加密因输出工作簿时的内存压力主动中断，退出码130；中断记录与原日志保留，没有完整轨迹，未将其列为通过。已完成网格以及同网格一维比较足以将秒级几何影响与数秒径向偏差区分。

## 文件与复算

积分通量主验证的实际源码为 `../source/geometry_solver_executed.py`，SHA256 为 `7f4fc1ecad068ff3f5af12f48374ec742b12bb5f85aa430d3c16c2cb8e24d80d`。`geometry_solver.py`仅增加冷启动目录接口；AST核对表明物性与离散算子完全未改。

先完成的Q1与调和通量粗算使用添加积分通量分支前的代码，原执行哈希为 `ae8887e12f10cf81375678ccdec66d155823f157673aa3f769c84a75713bc7eb`；当前程序的 `--flux harmonic`保留该离散。具体执行版本和磁盘哈希记录时点见 `source_execution_audit.json`。

所有对比数来自 `../validation/geometry/geometry_comparison.json` 与原始 NPZ。每份 NPZ 包含温湿二维快照、实际半径、参考径向坐标、轴向坐标、中截面/端面轨迹、全域最大值与位置以及事件前后完整状态。文件哈希见同目录清单。

从新目录复算的例子如下；若已有同名 JSON 或 NPZ，程序拒绝覆盖。输出接口冒烟试验仅验证路径与拒覆盖行为，不充当任何物理精度试验。

```bash
OPENBLAS_NUM_THREADS=1 python source/geometry_solver.py --case q23 --nr 80 --nz 128 --out-dir reproduced/q23_geometry_fresh --input-dir inputs
OPENBLAS_NUM_THREADS=1 python source/geometry_solver.py --case q23 --nr 80 --nz 0 --out-dir reproduced/q23_geometry_fresh --input-dir inputs
```
'''
(v1/'q123_closeout/一维二维适用范围.md').write_text(paragraph)

g=val('q4_80x128_integral' if 'q4_80x128_integral' in data else 'q4_40x128_integral')
q4text=f'''# 第四问：独立有限体积与二维验证

本轮只改变几何维数，不改变附录4物性、换热换质系数、空气等效边界或未来平台。长度固定为0.25m，材料均匀径向收缩，网格随材料运动。二维算例从初态完整联立求解温度与含水率，没有潜热项或等温投影。

## 参考区域算子

令 `ξ=r/R(t)`，材料速度与网格速度均为 `r R′/R`，相对速度为零。热方程和干基含水率方程分别为：

\\[
\\rho(C)c_p(C)\\partial_t T
=\\frac1{{R(t)^2\\xi}}\\partial_\\xi(\\xi k(C)\\partial_\\xi T)+\\partial_z(k(C)\\partial_zT),
\\]
\\[
\\partial_t C
=\\frac1{{R(t)^2\\xi}}\\partial_\\xi(\\xi D(T,C)\\partial_\\xi C)+\\partial_z(D(T,C)\\partial_zC).
\\]

以初始网格的对偶体积作权重，径向内部导通系数乘 `(R0/R)^2`，侧面交换乘 `R0/R`，轴向内部项与端面交换不变。这与当前体积、界面面积同时更新完全等价；没有重复加入收缩对流或干基C浓缩项。

## 实际完成的配对二维证据

{table(q4keys)}

默认积分配置在40×64→40×128时，配对几何差改变约 `0.0004711 s`；40×128→80×128改变约 `-0.0004190 s`。紧40×64配置采用 `rtol=2e-11、max_step=60s`，几何差相对默认配置变化约 `0.00564 s`，需纳入小差异的解释。

综合这些实际配对，端部对最湿点临界的影响约为提前 `0.0286–0.0343 s`，不宜把小数后所有位都解释为物理精度。最细80×128原场的全域最湿点为 `(r,z)=(0,0)`，并保存了等号根与前后1s的完整二维状态。

最细已完成网格全程中截面采样最大差为 `ΔT={g['max_mid_difference_T_C'][0]:.8g} K`、`ΔC={g['max_mid_difference_T_C'][1]:.8g} kg/kg`；端面差为 `{g['max_end_difference_T_C'][0]:.8g} K`、`{g['max_end_difference_T_C'][1]:.8g} kg/kg`。一维适用依据限定于中截面和最湿点。

## 径向离散误差与全域估计的边界

二维80×128原临界 `{data['q4_80x128_integral']['summary']['event_2d_s']:.9f} s` 晚于谱一维临界，是粗径向离散偏差；同网格一维也同样偏晚。不能宣称这个原始粗二维场在正式 `{q4['execution_s']:.2f} s` 时已经实际过线。

将各同网格几何差加到精细一维临界 `{q4['critical_s']:.9f} s` 上，得到全域临界估计 `{q4est[0]:.6f}–{q4est[1]:.6f} s`。结合一维临界斜率 `{slope:.12g} (kg/kg)/s`，正式执行时刻的全域最大C估计为 `{c4est[0]:.12f}–{c4est[1]:.12f}`。

上一段是“配对径向校正＋临界附近线性化”的数值估计，并非正式时刻的新二维场求值或严格上界。它与一维实测执行最大值 `{q4['execution_max_C']:.12f}` 一起支持一维工程执行时刻的适用性和保守方向；不宣称得到无限精度的二维最短时长。

独立积分FV的n320/n640紧配置临界为 `183931.398320243 / 183931.248715026 s`，二阶外推 `183931.198846621 s`。n640进一步减小步长与相对容差后为 `183931.259744118 s`；加上同一空间修正约为 `183931.209876 s`，与谱主临界相差约0.00184s。

上述空间外推和误差变化是实证数值检查，不是严格数学界。完整单网格数据与容差保存在各JSON中；没有把粗网格根直接替换为谱根，也没有为缩短时长调整物性、交换系数或观测收缩数据。

## 验证与来源

`jacobian_flux_checks.json`保存两种通量、三组物性的解析Jacobian方向差分检查；相对范数差约1e-10。`dimension_reduction_check.json`保存关闭端面交换、轴向均匀状态下二维与一维RHS的一致性；该结果只验证维数退化，不验证开放端部的几何误差。

`source_execution_audit.json`记录原始执行源与后续接口修订，`interface_revision_checks.json`记录AST不变及冷启动拒覆盖检查。任何旧JSON中的源码哈希都应结合该执行记录解释，因为旧实现曾在运行结束时才读取磁盘上的文件哈希。

原始矩阵、未舍入轨迹、完整二维事件场、配对对比JSON和CSV均在本目录。`corrected_global_event_estimates.json`清楚标识校正估计；两张PNG直接由保存的二维场生成，图中显示的事件场属于对应原始离散网格。
'''
(DEST/'README.md').write_text(q4text)
print('Wrote closeout geometry discussion and independent geometry README.')
