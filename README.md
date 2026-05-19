# RNA-seq 可重复分析论文复现作业

## 📌 项目基本信息
- **复现论文**：Designing RNA sequencing experiments: A practical guide to reproducible gene expression analysis
- **论文DOI**：10.1016/j.csbj.2025.12.015
- **论文链接**：https://doi.org/10.1016/j.csbj.2025.12.015
- **作业提交仓库**：https://github.com/D2RS-2026spring/RNA-seq-Reproduction
- **运行环境**：R 4.6.0 + Positron IDE

---

## 🧬 项目背景与目标
本项目是上述RNA-seq标准化分析指南论文的作业复现，核心目标为：
1.  掌握论文推荐的RNA-seq全流程分析逻辑
2.  验证从数据处理到差异分析、可视化的可重复性
3.  生成可直接用于作业提交的完整分析脚本与结果文件

---

## 📂 仓库结构说明
RNA-seq-Reproduction-1/
├── code/ # 核心分析脚本
│ └── run_all.R # 一键复现全流程的主脚本
├── count/ # 定量结果文件
│ └── counts.txt # 模拟 RNA-seq 计数矩阵
├── results/ # 分析结果与可视化图
│ ├── de_result.csv # 差异表达基因分析结果
│ ├── pca.png # 样本 PCA 聚类图
│ ├── volcano.png # 差异基因火山图
│ └── heatmap.png # 差异基因表达热图
└── README.md # 项目说明文档
