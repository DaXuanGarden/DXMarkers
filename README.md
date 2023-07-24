# DXMarkers：一键检索 "CellMarker"、"CellMarker2.0" 及 "PanglaoDB" 数据库的标记基因

在单细胞 RNA 序列分析 (scRNA-seq) 中，识别细胞类型是一项关键步骤。此过程常依赖于标记基因，即在某特定类型的细胞中高表达的基因。DXMarkers 正是以此为目标所设计的一款 R 包，它能够高效地在三个主流的标记基因数据库 "CellMarker"、"CellMarker2.0" 及 "PanglaoDB" 中检索特定基因，同时列出与之对应的细胞类型。此外，DXMarkers 还支持根据物种和组织类型对搜索结果进行筛选，从而提供更精确的信息。

## 安装指南

首先，你需要在 R 环境中安装 DXMarkers 包。以下是在本地环境中安装该包的操作流程：

```R
# 设定工作路径
setwd("/home/data/t050446/01 Single Cell Project/Acute Pancreatitis")

# 加载 devtools 包，用以安装本地 R 包
library(devtools)

# 安装本地版 DXMarkers 包
install_local("/home/data/t050446/Non-tumor research/DXMarkers_1.0.tar.gz")
```


## 一键检索标记基因

在进行 scRNA-seq 数据分析并鉴定出各个细胞群体的标记基因后，你可以使用 DXMarkers 的 `annotate_markers` 功能，在上述三个数据库中检索，同时获取每个基因对应的细胞类型：

```R
# 加载 DXMarkers 包
library(DXMarkers)

# 读取标记基因数据
top10_genes_data <- read.csv("top10_sorted_genes.csv", stringsAsFactors = FALSE)

# 利用 annotate_markers 功能进行基因注释，以小鼠胰腺为例
DXMarkers_result <- annotate_markers(top10_genes_data, "Mouse", "Pancreas")
```


## 查阅内建数据源

DXMarkers 亦提供了一项查阅内建数据源的功能，使你能够查看 "CellMarker"、"CellMarker2.0" 或 "PanglaoDB" 中的完整数据集：

```R
# 加载 DXMarkers 包
library(DXMarkers)

# 指定数据源
data_source <- "CellMarker"

# 浏览指定数据源的数据
data <- view_data_source(data_source = data_source)
```


DXMarkers 的开发目标是协助用户在大规模的单细胞 RNA-seq 数据中寻找关键的标记基因，同时能迅速且准确地检索这些标记基因在各数据库中的对应信息。无论你是单细胞数据分析的初学者，还是经验丰富的研究者，DXMarkers 都能帮助你节省时间，提升工作效率。