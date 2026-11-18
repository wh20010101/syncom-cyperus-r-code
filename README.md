# 论文配套 R 分析代码（Code-only）

本仓库仅包含论文分析所用的 R 代码，**不包含原始数据、中间数据和输出图/表**。
数据可用性见下文；运行前请按“输入数据（不包含在仓库中）”准备文件。

## 目录结构（仅代码）

```
.
├── README.md
├── LICENSE
├── .gitignore
├── 01_mantel/                 # 模块 1：Mantel 检验
│   └── 01_mantel.R
├── 02_metabolite_physiology/             # 模块 2
│   └── 02_metabolite_physiology_network.R
├── 03_microbe_metabolite/                         # 模块 3
│   └── 03_microbe_metabolite_network.R
├── 04_microbe_physiology/             # 模块 4
│   └── 04_microbe_physiology_network.R
├── 05_cooccurrence_rf/                     # 模块 5
│   ├── cooccurrence_network.R
│   └── random_forest.R
└── 06_factorial_anova/                               # 模块 6：双因素方差分析
    └── 06_factorial_anova.R
```

## 模块与运行顺序

| 顺序 | 模块 | 脚本 | 主要输入 |
|---|---|---|---|
| ① | 01_mantel（Mantel） | `01_mantel.R` | `SS/NSS` 的微生物、代谢物、生理指标 xlsx |
| ② | 代谢物—生理指标相关性网络 | `02_metabolite_physiology_network.R` | `metabolite_SS.csv / metabolite_NSS.csv` + 模块 1 的生理指标 xlsx |
| ③ | 微生物—代谢物相关性网络 | `03_microbe_metabolite_network.R` | `data/microbial_OTU_table.xls` + 模块 2 的代谢物 CSV |
| ④ | 微生物—生理指标相关性网络 | `04_microbe_physiology_network.R` | `data/microbial_OTU_table.xls` + `physiology.xlsx`、`groups.xlsx` |
| ⑤ | 微生物共现网络与随机森林 | `cooccurrence_network.R`、`random_forest.R` | `data/microbial_OTU_table.xls` |
| ⑥ | 06_factorial_anova | `06_factorial_anova.R` | `factorial_traits.xlsx` |

> 模块 5：4 个处理组（CK、S、SC1、SC5），SS 与 NSS 合并，每组 6 个样本；共现网络采用 `|r| > 0.7`、BH-FDR `P < 0.05`、前 100 属；随机森林采用前 30 属、`mtry = 6`、`ntree = 2000`、`set.seed(42)`。
>
> 模块 3 的微生物–代谢物网络取前 50 属；模块 4 的微生物–生理指标网络取前 100 属。

运行顺序：**① → ② → ③**；④、⑤、⑥ 在数据准备好后可独立运行。
模块 2、3 之间使用 `../` 相对路径读取中间表；模块 3、4、5 从 `../data/` 读取微生物丰度表。请保持目录层级不变。

## 输入数据（不包含在仓库中）

原始数据不随代码公开。若你已获得数据使用许可，请按下列位置放置文件：

| 文件 | 放置位置 | 说明 |
|---|---|---|
| 微生物丰度表 | `data/microbial_OTU_table.xls` | 制表符分隔；行 = OTU/种及分类注释，列 = 样本（`W*_*`）；脚本按 `genus` 聚合到属 |
| 代谢物筛选表 | `02_metabolite_physiology/metabolite_SS.csv`、`metabolite_NSS.csv` | 行 = 代谢物，列 = 注释信息 + 样本丰度（`W*_*`） |
| 生理指标（模块 1） | `01_mantel/SS_physiology.xlsx`、`NSS_physiology.xlsx` | 第一列 `Sample`，其余列 = 生理指标 |
| 微生物/代谢物（模块 1） | `01_mantel/SS_microbiome.xlsx`、`SS_metabolome.xlsx`、`NSS_microbiome.xlsx`、`NSS_metabolome.xlsx` | 第一列 `Sample`，其余列 = 特征 |
| 生理指标（模块 4） | `04_microbe_physiology/physiology.xlsx` | 第一列 `Sample`，其余列 = 生理指标 |
| 分组信息（模块 4） | `04_microbe_physiology/groups.xlsx` | 至少包含 `Sample`、`Group` 两列；请使用 `SS-CK` / `SS-S` / `SS-SC1` / `SS-SC5` 与对应 `NSS-...` 的标准写法 |
| 06_factorial_anova数据 | `06_factorial_anova/factorial_traits.xlsx` | `Sheet1`：42 行 × 39 列，Sample + 38 个指标 |

> 模块 6 使用全模型（42 个观测，误差 df = 28）：输出三效应 F/P/偏 η²/95% CI、Levene 与 Shapiro 检验，以及每个土壤背景内 7 个处理的 Tukey HSD（全模型 MSE，df = 28）和 Piepho 显著性字母；输出 `output/1_ANOVA_effects.xlsx` 与 `output/2_treatment_means_letters.xlsx`。原始数据不包含在仓库中。

## 运行方式

- **RStudio**：打开对应 `.R` 文件，点击 **Source**（或全选后 Ctrl+Enter）。
- **命令行**（在仓库根目录）：

```bash
Rscript "01_mantel/01_mantel.R"
```

大多数脚本开头有“自动切换到脚本所在目录”的路径适配代码；它只解决工作目录问题，**不代表数据已包含在仓库中**。数据缺失时脚本会报“文件不存在”，请先按上表放置数据。

## 环境依赖

R ≥ 4.1。脚本使用的包：

```
readxl  openxlsx  tidyverse  linkET  igraph  pheatmap
svglite  randomForest  writexl  rstudioapi
Hmisc  RColorBrewer
```

```r
install.packages(c(
  "readxl", "openxlsx", "tidyverse",
  "igraph", "pheatmap", "svglite",
  "randomForest", "writexl",
  "rstudioapi", "remotes",
  "Hmisc", "RColorBrewer"
))

# linkET 若不在 CRAN，可用 GitHub 安装
if (!requireNamespace("linkET", quietly = TRUE)) {
  remotes::install_github("Hy4m/linkET")
}
```

说明：

- 模块 1 使用 `readxl` / `linkET` / `tidyverse` / `openxlsx`；模块 6 使用 `readxl` + `writexl`，`car` 可选（平衡设计下 Type I/III 等价）。
- 建议在本地使用 `renv` 锁定包版本；本仓库暂未附带 `renv.lock`。
- 出图脚本使用 `Times New Roman`；Linux 环境如缺少该字体，请安装字体或把 `system_fonts` 改为可用的 serif 字体。

## 数据可用性

**English:**

> Raw sequencing data have been deposited in the NCBI Sequence Read Archive (SRA) under BioProject PRJNA1510944. Metabolomics data have been deposited in the China National Genomics Data Center (CNGDC) under BioProject PRJCA071506. Processed input tables used by the R scripts are not included in this repository and are available from the corresponding author upon reasonable request.

**中文：**

> 原始测序数据已存入 NCBI 序列读取档案（SRA）数据库，BioProject 编号为 PRJNA1510944；质谱代谢组学数据已存入中国国家基因组数据中心（CNGDC），BioProject 编号为 PRJCA071506。R 脚本使用的处理后输入表不包含在本仓库中，可在合理请求下向通讯作者获取。

## 代码可用性

> The R code used for the analyses is available at
> https://github.com/wh20010101/syncom-cyperus-r-code (DOI: [10.5281/zenodo.24353767](https://doi.org/10.5281/zenodo.24353767)).



## 许可

代码采用 MIT License，见 `LICENSE`。本仓库不包含数据，因此不附带数据许可；如果将来添加数据或元数据，请单独补充 CC BY 4.0 / CC0 等数据许可。


