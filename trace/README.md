古人类基因渗入检测分析流程：从 VCF 到 TRACE 精细定位
本项目基于 tsinfer、tsdate 和 TRACE 等开源工具，构建了一套从全基因组基因型数据（VCF）出发，检测古人类（尼安德特人、丹尼索瓦人）基因渗入片段的标准化、可复现的生物信息学分析流程。

特别说明：本项目不开发新算法，而是对已有工具的流程化整合与封装，旨在提高分析效率、可复现性和可读性。

项目概述
项目目标
搭建端到端的自动化流程（数据预处理 → ARG 推断 → 时间定标 → 渗入检测 → 可视化）。

以 1000 Genomes 第 22 号染色体 数据为例，验证流程稳定性，并撰写详细文档。

提供 Conda 环境文件（environment.yml）记录所有依赖，确保跨平台可复现。

提供结果可视化指南（群体 burden 柱状图、全基因组 Manhattan 图、IGV 交互式浏览）。

数据来源
数据集：1000 Genomes Project (Phase 3)

参考基因组：GRCh37 (hg19)

示例染色体：22 号染色体（可推广至其他染色体）

原始文件：chr22.vcf.gz（BGZF 压缩 VCF 格式）

样本数：2504 个样本（26 个人群，5 个超级群体：AFR, AMR, EAS, EUR, SAS）

数据来源：1000 Genomes FTP 或本地路径

用户可根据需要替换为其他 VCF 文件（需符合相同格式要求）。

环境与依赖
操作系统与基础环境
工具	说明	安装方式
Linux / WSL2	所有分析均在 Linux 环境下运行（推荐 Ubuntu 20.04+）	自行安装（Windows 用户可通过 WSL2 安装）
Bash	脚本运行环境（#!/usr/bin/env bash）	系统自带（Linux/WSL）
Conda	环境管理与包依赖管理（≥23.0.0）	Miniconda 安装指南
R (≥4.0)	用于部分可视化（如 ibdmix.R）	conda install -c conda-forge r-base
Python (3.10.x)	脚本运行环境（tsinfer、TRACE 等均基于 Python）	conda create -n gu python=3.10
核心软件与 Python 库
工具名称	版本建议	功能描述	安装方式
bcftools	≥1.16	VCF 文件处理（拆分、排序、索引、过滤）	conda install -c bioconda bcftools
vcf2zarr (bio2zarr)	≥0.2.0	将 VCF 转换为 Zarr 格式（.vcz），供 tsinfer 高效读取	pip install bio2zarr
tsinfer	≥0.6.0	从基因型数据推断祖先重组图（ARG）	pip install tsinfer
tsdate	≥0.3.0	对 ARG 进行时间定标，生成带时间的树序列	pip install tsdate
tskit	≥0.5.0	处理树序列文件（加载、简化、保存）	pip install tskit
TRACE	latest	从定标 ARG 中检测古人类基因渗入片段	pip install git+https://github.com/YulinZhang9806/trace.git
zarr	≥2.13.0	处理 Zarr 格式数据（添加 REF_allele 等）	pip install zarr
numpy	≥1.23.0	Python 数值计算（辅助脚本）	pip install numpy
pandas	≥1.5.0	数据处理与表格操作	pip install pandas
matplotlib	≥3.5.0	生成可视化图表（折线图、柱状图等）	pip install matplotlib
plotly（可选）	≥5.0.0	交互式 ARG 浏览器（Lorax）的支持库	pip install plotly
所有依赖均可通过提供的 environment.yml 文件一键创建 Conda 环境（conda env create -f environment.yml）。

完整分析流程
步骤 1：数据准备
下载 1000 Genomes 的 VCF 文件（chr22.vcf.gz）及索引（.tbi）。

准备样本群体信息文件（samples_v3.ALL.panel），包含三列：sample、pop、super_pop。

可使用项目提供的 gen.clean.sh 脚本自动完成。

步骤 2：构建祖先重组图（ARG）
使用 tsinfer 和 tsdate 从 VCF 推断时间定标的 ARG：

VCF → Zarr：

bash
vcf2zarr convert chr22.vcf.gz chr22_full.vcz
添加 REF_allele 数组（用于指定祖先状态）
（可通过脚本自动完成）

运行 tsinfer.infer() 推断初始 ARG。

运行 tsdate.date() 进行时间定标，生成 .trees 文件。

关键点：确保 ARG 的个体元数据中包含样本 ID（tsinfer 默认从 VCF 继承），或额外生成 sample_map.tsv 供 TRACE 使用。

步骤 3：TRACE 精细分析
TRACE 工作流分为四个子步骤：

子命令	功能
trace-extract	从 ARG 中提取每个节点在指定区域的祖先“观察值”
trace-infer	推断每个节点的祖先状态（后验概率）
trace-summarize	汇总连续的渗入片段（基于后验和长度阈值）
trace.py combine	合并所有节点片段，生成窗口 burden 表（trace_all_gw.csv）和样本长度表（trace_all_length.csv）
精细定位设置：

目标区域：chr22:30,000,001–31,000,000

窗口大小：50 kb

使用 TRACE_INCLUDE_DIR 传入 BED 文件限定区域，TRACE_WINDOW_SIZE=50000 提高分辨率。
