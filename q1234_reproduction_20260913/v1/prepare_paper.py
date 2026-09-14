"""Prepare an updated appendix without changing the accepted paper body."""
from pathlib import Path
import shutil

BASE = Path(__file__).resolve().parent
REPO = BASE.parents[1]
OLD = REPO / "paper/q1234_draft_v1"
SOURCE = BASE / "paper_source"
PACKAGE = BASE / "package/支撑材料"


def main():
    if SOURCE.exists():
        raise FileExistsError(SOURCE)
    shutil.copytree(PACKAGE, SOURCE)
    for filename in ["main.tex", "ai_details.tex", "validate_paper.py"]:
        shutil.copyfile(OLD / filename, SOURCE / filename)
    shutil.copytree(OLD / "sections", SOURCE / "sections")
    shutil.copytree(OLD / "tables", SOURCE / "tables")
    (SOURCE / "evidence").mkdir()
    shutil.copyfile(OLD / "evidence/table_expected.json", SOURCE / "evidence/table_expected.json")
    appendix = SOURCE / "sections/appendix.tex"
    text = appendix.read_text(encoding="utf-8")
    first_code = text.index(r"\sourcecode{匿名快照启动与冒烟检查}")
    introduction = r"""\section{支撑材料与完整采用代码}
\label{sec:complete-code}
支撑包包含完整求解源码、题面附件、提交结果与数值复现入口。下列原求解代码保留；新增的完整复现控制、常规Excel导出、逐格核对及数值图生成代码一并列于本附录。代码附录不计入正文页数。

\codegroup{完整复现入口与验收范围}
在新的Python环境中安装支撑包的\texttt{requirements.txt}，运行\texttt{python -X utf8 reproduce.py --out ../reproduced --jobs 2}。程序仅以包内源码和题目附件重新计算四问主解、数值与二维对照、交换系数及后段环境情景，再核对论文表格、终点和提交结果。包内已有结果仅供比较，不作为求解初态。

输出目录包含新算的未舍入场、实际命令及退出码、独立工作簿读回、比较记录和数值核查图。第二问新生成的\texttt{result2.xlsx}包含完整全过程；压缩提交包中的同名小工作簿保存前三小时，全过程另以\texttt{NPZ}提供。新导出入口使用Python标准库与openpyxl，无需专有工作簿运行时。旧版导出代码保留以供追溯。

六张题设表按论文显示精度逐格比较；临界根和未舍入量使用记录在核验程序中的数值容差。复现结果不改变模型的条件性，也不将第四问一维执行值提升为有限圆柱连续全域的严格保证。运行方法及结果见支撑包的复现说明与验收记录。

\codegroup{原求解与诊断源码}
"""
    new_code = "\n\\codegroup{新增完整复现与常规导出代码}\n"
    paths = ["reproduce.py", *[p.relative_to(PACKAGE).as_posix() for p in sorted((PACKAGE / "reproduction").glob("*.py")) if p.name != "__init__.py"]]
    for path in paths:
        title = Path(path).stem.replace("_", r"\_")
        new_code += r"\sourcecode{" + title + "}{" + path + "}\n\n"
    appendix.write_text(introduction + text[first_code:] + new_code, encoding="utf-8")
    ai = SOURCE / "ai_details.tex"
    original = ai.read_text(encoding="utf-8")
    supplement = r"""
\section{支撑包复现链补充记录}
本轮根据“支撑材料中的代码能否复现论文结果”及“你来完成”的请求，补充独立工作目录中的实际数值复跑、原始输入与源码哈希核对、常规Excel导出、逐格比较以及运行说明。保留原模型和数值算子，新增复现控制与核对代码。

原题设表、临界根和完整提交数据分别按显示精度与记录的数值容差核对。程序读回导出工作簿，并以故意改错的样本检查比较器能够报警。实际结果、运行命令与退出状态存于复现验收记录。

以上为工具实际执行的技术核查；新增代码及记录的最终采纳和人工审阅由参赛队负责，自动检查不替代人工确认。原文模型适用范围和第四问二维全域验证限制保留。
"""
    ai.write_text(original.replace(r"\end{document}", supplement + "\n" + r"\end{document}"), encoding="utf-8")
    print(SOURCE)


if __name__ == "__main__":
    main()
