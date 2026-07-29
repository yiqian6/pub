1\. 项目概述

1.1 项目目标

本项目旨在整合现有开源工具，构建一套从 全基因组基因型数据（VCF） 出发，最终检测古人类（如尼安德特人、丹尼索瓦人）基因渗入片段的标准化、可复现的生物信息学分析流程。



具体工作包括：



搭建端到端的自动化流程：将数据预处理、ARG 推断、时间定标、渗入检测、可视化等步骤串联成一个完整的、可一键运行的工作流（Shell + Python 脚本），降低用户手动干预和出错的风险。



验证与文档化：以 1000 Genomes 第 22 号染色体数据为例，测试整个流程的稳定性和可复现性，并编写详细的项目文档（README），包括每一步的输入、输出、命令解释和预期结果，方便他人理解和使用。



提供可复现的环境配置：通过 Conda 环境文件（environment.yml）记录所有软件依赖和版本，确保流程在不同计算环境下都能顺利运行。



结果可视化指南：提供群体承载率柱状图、全基因组 Manhattan 图以及交互式 ARG 浏览器（Lorax）的使用示例，帮助用户直观解读分析结果。



特别说明：本项目不开发新的算法或软件，而是对 tsinfer、tsdate、TRACE、bcftools、vcf2zarr 等已有开源工具的流程化整合与封装，旨在提高分析效率、可复现性和可读性。



1.2 数据来源

本项目使用 1000 Genomes Project (Phase 3) 的公开数据作为测试案例。



参考基因组版本：GRCh37 (hg19)



染色体：第22号染色体（chr22）作为示例，流程可推广至其他染色体。



原始文件：chr22.vcf.gz（BGZF 压缩的 VCF 格式）



样本数：2504 个样本（来自 26 个人群，分为 5 个超级群体：AFR, AMR, EAS, EUR, SAS）



数据来源：官方 FTP 站点或通过 bcftools 从本地路径读取。



用户可根据需要替换为其他 VCF 文件（需符合相同格式要求）。



1.3 软件工具与依赖

本项目依赖以下主要开源软件和 Python 库，所有工具均可在 Linux 环境（如 WSL2）下运行。



工具名称	版本建议	功能描述	安装方式

bcftools	≥1.16	VCF 文件处理（拆分、排序、索引、过滤）	conda install -c bioconda bcftools

vcf2zarr (bio2zarr)	≥0.2.0	将 VCF 转换为 Zarr 格式（.vcz），供 tsinfer 高效读取	pip install bio2zarr

tsinfer	≥0.6.0	从基因型数据推断祖先重组图（ARG）	pip install tsinfer

tsdate	≥0.3.0	对 ARG 进行时间定标，生成带时间的树序列	pip install tsdate

tskit	≥0.5.0	处理树序列文件（加载、简化、保存）	pip install tskit

TRACE	latest	从定标 ARG 中检测古人类基因渗入片段	pip install git+https://github.com/YulinZhang9806/trace.git

zarr	≥2.13.0	处理 Zarr 格式数据（添加 REF\_allele）	pip install zarr

numpy	≥1.23.0	Python 数值计算（辅助脚本）	pip install numpy

python	3.10.x	脚本运行环境	conda create -n tsinfer\_env python=3.10

conda	≥23.0.0	环境管理与包依赖管理	自行安装

