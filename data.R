######数据读取并保存#######
# 设置工作目录
setwd("/home/data/t050446/Non-tumor research/DXMarkers")

# 定义一个函数，用于读取CSV文件并保存为RDS格式
read_and_save_as_rds <- function(csv_path, rds_filename) {
  data <- read.csv(csv_path)
  # 保留PMID列，删除其右侧的所有列
  data <- data[, 1:which(names(data) == "PMID")]
  saveRDS(data, file = rds_filename)
}

# 读取并保存第一个文件
cellmarker1_path <- "data/CellMarker_all_cell_markers.csv"
cellmarker1_rds <- "data/CellMarker_all_cell_markers.Rds"
read_and_save_as_rds(cellmarker1_path, cellmarker1_rds)

# 读取并保存第二个文件
cellmarker2_path <- "data/CellMarker2.0_Cell_marker_All.csv"
cellmarker2_rds <- "data/CellMarker2_0_Cell_marker_All.Rds"
read_and_save_as_rds(cellmarker2_path, cellmarker2_rds)

# 读取并保存第三个文件
cellmarker3_path <- "data/PanglaoDB_markers_27_Mar_2020.csv"
cellmarker3_rds <- "data/PanglaoDB_markers.Rds"
read_and_save_as_rds(cellmarker3_path, cellmarker3_rds)
# 加载RDS格式的数据集
CellMarker_all_cell_markers <- readRDS("data/CellMarker_all_cell_markers.Rds")
CellMarker2_0_Cell_marker_All <- readRDS("data/CellMarker2_0_Cell_marker_All.Rds")
PanglaoDB_markers <- readRDS("data/PanglaoDB_markers.Rds")

######数据合并######

# 设置工作目录
setwd("/home/data/t050446/Non-tumor research/DXMarkers")

# 加载已经保存的RDS格式的数据集
CellMarker_all_cell_markers <- readRDS("data/CellMarker_all_cell_markers.Rds")
CellMarker2_0_Cell_marker_All <- readRDS("data/CellMarker2_0_Cell_marker_All.Rds")
PanglaoDB_markers <- readRDS("data/PanglaoDB_markers.Rds")

# 为每个数据集添加一个来源标识列
CellMarker_all_cell_markers$source <- "CellMarker"
CellMarker2_0_Cell_marker_All$source <- "CellMarker2.0"
PanglaoDB_markers$source <- "PanglaoDB"

# 将这三个数据集按列拼接
DXMarkers_data <- rbind(CellMarker_all_cell_markers, CellMarker2_0_Cell_marker_All, PanglaoDB_markers)

# 保存合并后的数据
saveRDS(DXMarkers_data, file = "data/DXMarkers_data.Rds")

#####分隔并保存#####
library(tidyverse)

# 设置工作目录
setwd("/home/data/t050446/Non-tumor research/DXMarkers")

# 加载已经保存的RDS格式的数据集
CellMarker_all_cell_markers <- readRDS("data/CellMarker_all_cell_markers.Rds")
CellMarker2_0_Cell_marker_All <- readRDS("data/CellMarker2_0_Cell_marker_All.Rds")
PanglaoDB_markers <- readRDS("data/PanglaoDB_markers.Rds")

# 提前对内置数据进行修改，将逗号分隔的基因名分开并复制其他列的信息
CellMarker_all_cell_markers <- CellMarker_all_cell_markers %>%
  separate_rows(geneSymbol, sep = ", ")

CellMarker2_0_Cell_marker_All <- CellMarker2_0_Cell_marker_All %>%
  separate_rows(geneSymbol, sep = ", ")

PanglaoDB_markers <- PanglaoDB_markers %>%
  separate_rows(geneSymbol, sep = ", ")

# 为每个数据集添加一个来源标识列
CellMarker_all_cell_markers$source <- "CellMarker"
CellMarker2_0_Cell_marker_All$source <- "CellMarker2.0"
PanglaoDB_markers$source <- "PanglaoDB"

# 将这三个数据集按列拼接
DXMarkers_data <- rbind(CellMarker_all_cell_markers, CellMarker2_0_Cell_marker_All, PanglaoDB_markers)

# 保存合并后的数据
saveRDS(DXMarkers_data, file = "data/DXMarkers_data.Rds")
