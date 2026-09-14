# 匿名可复现代码支撑

此目录保存正文采用的 42 个完整源码文件，另有启动脚本 `launch.py`、小型题面输入和必要模板。`project/` 镜像原项目相对路径。完整源码均输入论文附录；没有用节选代替核心算法。这里是代码与小输入支撑，较大的冻结结果数组、正式工作簿仍由独立结果交付提供。

本次制稿没有重跑四问全部 PDE，也没有重新导出全部工作簿。实际完成的检查见 `smoke_result.json`：64 个快照文件哈希、42 个 Python 文件 AST（含启动脚本）、25 个安全 `--help` 入口以及 1 个 JavaScript `node --check` 均通过。安全入口检查同时验证其正常 import 依赖；它不是独立数值求解、全量运行或物理校准。

## 1. 目录与运行环境

| 范围 | 核心路径 | 作用 |
|---|---|---|
| Q1 | `project/q1_complete_delivery/q1_delivery/` | 主求解、验证、加密、舍入、Excel、历史报告链 |
| Q2 机制与情景 | `project/q2_refinement_delivery/` | 对已存场分解扩散系数作用，并保留生成同网格 320 FV 参数情景的实际代码 |
| Q2/Q3 现行统一轨迹 | `project/q4_complete_delivery/v1/q123_closeout/source/` | 全程每秒轨迹及与旧前三小时、Q3 的核对 |
| Q3 | `project/q3_refinement_delivery/source/` | 严格执行点和 60 s 轨迹、工作簿导出与验证 |
| Q4 与二维 | `project/q4_complete_delivery/v1/source/` | 收缩主模型、数值预算、执行点、有限圆柱及工作簿 |
| 本轮复审 | `project/q1234_overall_review_delivery/v1/` | 密度/能量代数诊断、4 h 环境续算、二维局部终点与独立证据核验 |
| 本稿制图 | `project/paper/q1234_draft_v1/build_assets.py` | 6 张表、5 幅图的实际只读后处理 |

本轮 smoke 实际 Python 为 3.13.14，NumPy 2.3.5、SciPy 1.17.0、Matplotlib 3.11.0、openpyxl 3.1.5；无需安装专有工作簿工具即可通过这些检查。上述版本是本轮入口检查环境，不表示全部历史数值结果都由这一环境生成。各原始子目录 requirements 保留供核对；它们的绘图库固定版本并不完全相同，应分开建立环境，避免将一次冒烟当作跨版本逐位复现证明。

可在新虚拟环境安装 `requirements-smoke.txt` 后，从论文目录运行：

```console
python support/launch.py smoke --log support/smoke_result.json
python support/launch.py list
```

若 Node 已可用，可加 `--node node` 进行 JavaScript 语法检查。`node --check` 不加载 `@oai/artifact-tool`，因此不证明该专有包已安装或导出流程能完成。Q1/Q3 提供标准库 OOXML 导出路径；Q4 的分块导出需要 `@oai/artifact-tool`、环境变量 `CODEX_PRIMARY_RUNTIME_NODE` 和 `CODEX_PRIMARY_RUNTIME_NODE_MODULES`。可用 `DRYING_WORKBOOK_TEMP` 指定匿名临时目录；目录符号链接权限、Linux 内存/字体探针、字体渲染均需在目标环境另验。本轮没有执行此导出链。

## 2. 明确启动的数值计算

下面的 `run` 命令会积分 PDE，耗时与内存可显著增加。本轮没有运行这些命令。启动脚本默认仅 `smoke`，数值命令必须显式给出一个尚不存在且位于 `support/project/` 以外的新输出目录。

```console
python support/launch.py run q1-main --out reproduced/q1
python support/launch.py run q23-main --out reproduced/q23
python support/launch.py run q3-main --out reproduced/q3
python support/launch.py run q4-main --out reproduced/q4
python support/launch.py run q2-sensitivity --out reproduced/q2_sensitivity
```

Q1 完整求解、验证、加密、工作簿和历史报告链可用 `run q1-full --out reproduced/q1_complete` 显式执行。Q1 报告链使用随包 `paper_template.md`，其文字是原阶段模板，不能替代本次四问复审后的正文。

| 任务 | 采用配置/输出 | 后续验收 |
|---|---|---|
| Q1 主解 | 10240 FV，题面输入，主代码默认严格容差及边界通量求积 | `q1_validate.py`、`q1_refine.py`、`q1_rounding_audit.py`；再 Excel 全单元格核验 |
| Q23 全程 | 80 节点 Radau，rtol=2e-13，max_step=30 s，输出步长 1 s，206906.76 s | `audit_unified_q23.py` 需参考谱解、历史 Q2 数组/工作簿和现行 Q3 |
| Q3 | 80 节点 Radau，rtol=2e-11，atol_C=2e-13，max_step=120 s；0.03 s 预算、0.36 s 向上量化 | 独立节点/方法/有限体积差与全表舍入；预算本身仍是经验数值证据 |
| Q4 | 120 节点 Radau，rtol=2e-11，max_step=90 s，线性半径与质量积分诊断 | 须先完成下面的 Q4 对照矩阵，才能 `finalize_q4.py` |
| Q2 参数情景 | `q2_tests.run_sensitivity`，同为 320 FV，rtol=1e-10，max_step=20 s | 7 个情景与同网格基线比较；hm/hT ±20% 不是置信区间 |

第三问现行工作簿来自冻结的 80 节点 Radau 轨迹。新入口采用同一算子和自动执行点；重新生成的微小未舍入差应按数值容差和四位单元格复核，不能要求新求解的压缩文件 SHA 与历史 NPZ 相同。

## 3. Q4 必要重建顺序与输出接口

先在新的工作副本中准备 `output/`、`validation/`、`validation/geometry/`。将主解命名为 `output/main.npz` 和 `main.json`。其余已采用对照配置为：

| 目标名 | 程序与关键参数 |
|---|---|
| `validation/spectral80` | `run_q4.py --n 80 --method Radau --rtol 2e-10 --max-step 180` |
| `validation/time120_bdf` | `run_q4.py --n 120 --method BDF --rtol 1e-11 --max-step 45` |
| `q4_320x0_integral` | `geometry_solver.py --case q4 --nr 320 --nz 0 --flux integral --rtol 2e-9 --max-step 240` |
| `q4_320x0_integral_rtol2e-10` | 同上，rtol=2e-10，max_step=120 |
| `q4_640x0_integral_rtol2e-10` | nr=640，nz=0，rtol=2e-10，max_step=120 |

谱程序需显式提供 `--out` 文件前缀。几何接口需 `--input-dir` 指向随包 `attachment1.xlsx`/`attachment2.xlsx` 所在目录，并为每种容差给不同的 `--out-dir`。几何程序默认标签不含容差，汇总时应在新工作目录按表中标签复制并保留原/新文件哈希，不能覆盖另一容差结果。

二维检查仍需要同一 nr 下 nz=0 与 nz>0 配对；`geometry_analysis.py` 读取全部配对实际输出，`geometry_report.py` 还读取 Q2 收尾汇总、Q4 主解和执行点。它们保留历史生成逻辑，其中关于“上界”等文字不能自动充当本次论文结论，必须采用正文已修订的“离散证据/校正估计/严格数学界”区分。`geometry_solver_executed.py` 是原实际求解算子快照；提供输出路径接口的 `geometry_solver.py` 用于新运行，两份均完整保留。

`check_numerics.py --base 工作副本` 读取主解和对照，并会写比较 JSON。之后 `finalize_q4.py --base 工作副本` 才从实际差异形成预算、选择执行点并产生表6/工作簿数据。它会确认执行点已在完整场中实算且 `max(C)<0.15`，若新求解的量化点未被存储，会拒绝导出；应再计算该点而不是补造一行。`package_workbooks.py --help` 给出 Q2/Q4 导出接口；只向新 `--out-dir` 写文件。

半径程序在 72 h 以后保持最后观测半径。正式 Q4 执行点在附件2时间范围内；超过 72 h 的新情景须先明确新的半径外推假设，不得把该默认延伸作为已观测输入。

## 4. 冻结数组后处理的外部材料

以下大文件未重复纳入本代码目录。需要复算历史图表/诊断时，应从独立结果交付核对哈希后按相同路径放入**工作副本**，或者通过程序提供的显式参数指向只读原证据。不能将缺少外部材料的脚本失败称为模型失败。

| 用途 | 必需外部路径/内容 |
|---|---|
| Q1 图表和密度诊断 | `q1_complete_delivery/q1_delivery/output/q1_unrounded.npz` 与 `table_temperature.csv`、`table_moisture.csv` |
| 当前 Q2/Q3 全程 | `q4_complete_delivery/v1/q123_closeout/output/q23_unified.npz` |
| Q3 表5及工作簿 | `q3_refinement_delivery/output/solution.npz`、`table5.csv`、`end_event.json`；`verify_frozen.py` 另需该交付的完整冻结清单/结果 |
| Q4 表6与收缩图 | `q4_complete_delivery/v1/output/main.npz`、`main.json`、`table6_unrounded.csv`、`end_event.json` |
| Q2 历史前三小时对照 | `q2_final_delivery/output/q2_unrounded.npz`、历史 `result2.xlsx`、两张 table CSV；此历史工作簿只作对照 |
| Q2 机制 | `q2_refinement_delivery/reused/q2_final_delivery/output/q2_unrounded.npz`、`validation/validation.json`、`reused/sensitivity_tracks.npz` |
| Q2 投影重建 | `prepare_reused.py` 的 `--handoff`（含原清单）和 `--arrays-zip`（7 个原参数数组）；由清单 SHA 核验后只选择 r=0/1/2 cm，无 PDE |
| 4 h 环境续算 | 上述 Q23 和 Q4 原始完整场；运行 `environment_probe.py --case q3/q4 --n 40 --scenario all --out 新目录` 会积分，非默认行为 |
| 独立环境核查 | `q1234_overall_review_delivery/v1/evidence/environment/` 的 14 对 JSON/NPZ，以及元数据引用的原始初态和原算子源码 |
| 二维局部终点 | `q4_complete_delivery/v1/validation/geometry/` 的 `q23_80x256_integral`、`q4_80x128_integral` JSON/NPZ；`audit_geometry_endpoints.py` 会局部续算 |
| 二维范围只读核查 | 同一 geometry 目录下脚本列出的各配对 JSON/NPZ；`audit_saved_geometry_scope.py` 顶层即执行，不能用 `--help` 试探 |
| 本稿全部表图 | 上述数据、`paper/q1_q2_stage/figures/q1_profiles`/`q2_profiles` 的 PDF/PNG、环境情景 JSON；运行镜像的 `build_assets.py`，生成其镜像论文目录内表图 |

Q2 `prepare_reused.py` 严格验证的是**原封存包**。若采用新运行的 320 FV 结果，应在新目录记录新的运行配置和源 SHA，并明确这是重算证据；不能伪装成原包哈希。`q2_mechanism.py` 可通过 `--raw --sensitivity --environment --validation --out` 显式接入已经确认的数据。

Q3 `export_workbook.py` 的 `--root` 需包含 `output/solution.npz` 与 `inputs/result3_blank.xlsx`。如使用 `launch.py run q3-main` 输出，应在新的工作目录整理成该结构，再调用 `--engine portable` 或 `--engine artifact`，随后运行 `validate_xlsx.py` 全表验证。

## 5. 匿名化、完整性和未验证边界

`snapshot_manifest.json` 同时记录原文件 SHA256、副本 SHA256、字节是否一致和文件作用，不含身份与提交标记。除少量来源说明、Git 元数据读取和工作簿临时目录参数化外，副本保留原字节；未修改方程、物性、时间节点、离散通量、收缩映射和阈值判据。内部审核另存逐行替换说明，不作为匿名投稿材料。

少量原源码含固定来源标记，匿名化必然改变其 SHA。旧环境 JSON 精确保存原算子 SHA，独立核查保留这一严格门槛；在匿名副本上核验旧 JSON 可能按设计失败。正确做法是先在内部原字节快照核验旧结果，或用匿名源码生成新情景与新哈希。**不要关闭断言，不要把原 SHA 填成匿名文件的实际 SHA。** 数值语义不变与文件字节一致是两个不同命题。

没有默认调用 `geometry_report.py`、`audit_saved_geometry_scope.py` 等顶层执行脚本；只对它们做 AST 检查。没有执行全套数值流程、环境 PDE、专有工作簿导出、冻结结果全量重审或原历史报告再生成。完整代码的保留使这些工作可以检查和有条件重建，不代表它们在本次编稿时已全部运行。


## 当前正文图表与提交包导览

[图表来源清单](图表来源清单.md)对应现行12幅正文图：四张分问流程图、一张问题一中截面与边界示意图、七张数值图。最新图内字号已按图3统一；当前成品图在提交包的 `../figures/` 目录，来源清单保留数组路径、哈希及各轮图形解释。

提交ZIP的顶层目录为 `支撑材料/`。`results/` 中有 `result1.xlsx` 至 `result4.xlsx`：其中 `result2.xlsx` 保存前3小时网格，完整全过程逐秒结果另存于同目录的 `result2_全程轨迹_206907x21.npz`。该NPZ包含 `time_s`、`radius_cm`、`temperature_TC`、`moisture_C`；温度与含水率数组均为206907×21。代码、输入和模板在本 `support/` 目录，图形在 `figures/`，AI说明和成员清单在顶层，参考文献在 `literature/`。论文PDF不包含在此支撑ZIP中。

图表来源清单中的 `../evidence/` 链接指向完整版LaTeX源稿的历史核验记录；精简提交ZIP未重复包含该目录。历史数据与数值内核未改，本文档所述历史检查不代表本次打包重新运行了求解或验证。

第四问二维终点的详细讨论已从正文移除，现有核查记录随包保留于 [二维终点记录](evidence/q4_endpoint_20260912/README.md)。其中JSON与NPZ均为已有结果的原字节副本，配有来源与哈希清单；正文第11.3节保留粗网格正式点略高于阈值、全域达标尚待进一步验证的限制说明。此次整理没有重新积分。
