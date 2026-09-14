# 网页模型阅读入口

本目录将项目中的4篇论文和原题提取为普通 Git 存储的 UTF-8 Markdown。
通过 GitHub 连接阅读时，先打开本文，再按链接读取目标论文全文或单页。
单页文件用于按页引用、分批阅读及避免长文件截断；它们不是摘要。

第一问已有新的 [独立审核入口](../review_q1_20260910/README.md) 和 [完成情况与答案审核报告](../review_q1_20260910/第一问审核报告.md)，其中链接本机重跑证据及独立热参考解。审查当前模型时请一并阅读，区分模型内数值验证与物理假设是否成立。

| 文档 | PDF页数 | 连续文本 |
|---|---:|---|
| [Estimation of thermo-physical properties of products with cylindrical shape during drying: The coupling between mass and heat](papers/da-silva-2014/README.md) | 9 | [全文](papers/da-silva-2014/fulltext.md) |
| [Combined Heat and Mass Transfer Associated with Kinetics Models for Analyzing Convective Stepwise Drying of Carrot Cubes](papers/chupawa-2022/README.md) | 19 | [全文](papers/chupawa-2022/fulltext.md) |
| [A Non-Isothermal Moving-Boundary Model for Continuous and Intermittent Drying of Pears](papers/adrover-2020/README.md) | 22 | [全文](papers/adrover-2020/fulltext.md) |
| [基于动网格的白萝卜热风干燥热质传递研究](papers/wu-2022/README.md) | 10 | [全文](papers/wu-2022/fulltext.md) |
| [2026年A题：药材的烘干问题](problem/README.md) | 4 | [全文](problem/fulltext.md) |

## 阅读与引用

- 如全文返回不完整，按每篇 README 的页码目录逐页读取；不要将一次截断响应视为已读全文。
- 页码是PDF物理页码。引用时注明论文、PDF页码和实际找到的章节/公式编号。
- `layout.txt`保留部分原始版面关系，可与阅读顺序版互相对照。
- 双栏的da Silva和中文白萝卜论文使用经抽查更连贯的`-raw`顺序，其余文档使用默认顺序；这不保证每张跨栏表的顺序完整。
- `⟦U+xxxx⟧`显式标出未可靠解码的字符。部分标记在原PDF中是负号、乘号或括号，不能删除后直接代入计算。
- 阅读前请查看[文本质量与公式风险说明](QUALITY.md)，其中列出了实际抽查范围和重点回看页。
- 数学公式的上下标、分式、希腊字母、双栏表格以及图像不能仅凭文本保证准确；关键公式回看原件或上传相关页面图像。
- 原始PDF全部保留。一篇约120 MiB的PDF经Git LFS保存；本目录的所有文本均为普通Git文件，不依赖LFS下载。
- 本目录补充的是当前全文访问证据。旧分析中“未获取全文”的记录反映当时的阅读状态，不代表现在仍缺原件；全文提取也不等于所有内容已由人逐式核验。

## 权限与文件格式分开检查

1. 先通过GitHub连接读取仓库根目录的README.md。普通文本也返回404时，先核对登录账号和仓库新名称；如果仓库是私密的，再核对连接对该仓库的授权。
2. 确认连接能读README后，再读取本目录的Markdown；不要要求文本读取接口把PDF二进制或LFS指针当成论文正文。
3. 仓库名是 `Hel10o/CUMCM26`（此前为 `Hel10o/libai`，旧地址由 GitHub 重定向），当前按所有者要求设为公开。如以后改为私密且授权范围按仓库选择，请在当前GitHub连接的配置中确认它已被包含，再重试读取。

## 再生成与完整性

安装Python包 `pypdf` 并确保Poppler的 `pdftotext` 在PATH中，在仓库根目录运行：

```bash
python scripts/export_readable_pdfs.py
```

`manifest.json`记录原件SHA-256、提取模式与页数、逐页字符数、低文本页、替换字符及控制字符标记数。
脚本核对原件哈希在提取前后不变；这些指标验证提取覆盖和文件完整性，不证明公式语义无误。
