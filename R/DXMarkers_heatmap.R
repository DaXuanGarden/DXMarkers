#' Generate Heatmap for DX Markers with log(*+1) transformed expression
#'
#' This function generates a heatmap for DX Markers from single-cell RNA sequencing data (scRNA) with log(*+1) transformed expression values, and saves it as a PNG file.
#'
#' @param scRNA An object containing single-cell RNA sequencing data. It can be Seurat, SingleCellExperiment, or other supported formats.
#' @param DXMarkers_result A list containing gene signatures for different cell types.
#' @param output_file Path to the output PNG file where the heatmap will be saved.
#' @param width Width of the output PNG image in pixels. Default is 800 pixels.
#' @param height Height of the output PNG image in pixels. Default is 800 pixels.
#' @param res Resolution (dots per inch) for the output PNG image. Default is 120 dpi.
#'
#' @name DXMarkers_heatmap
#' @import Seurat
#' @import pheatmap
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @export
#' @examples
#'
#' # 载入所需的包
#' library(Seurat)
#'
#' # 加载单细胞RNA数据对象
#' # Replace "scRNA" with the actual object containing your scRNA data
#' scRNA <- your_scRNA_data
#'
#' # 定义gene signatures的列表
#' # Replace "DXMarkers_result" with the actual list of gene signatures
#' DXMarkers_result <- your_DXMarkers_result
#'
#' # 指定输出的PNG文件名
#' output_file <- "heatmap_output.png"
#'
#' # 生成热图
#' DXMarkers_heatmap(scRNA, DXMarkers_result, output_file)
# 修改后的 DXMarkers_heatmap 函数
DXMarkers_heatmap <- function(scRNA, DXMarkers_result, output_file, width = 800, height = 800, res = 120) {

  # 调用 SCINA_markers 函数
  signatures <- SCINA_markers(DXMarkers_result)

  # 定义将列表转换为数据框的函数
  transform_signature_to_df <- function(signatures) {
    data_frame <- stack(signatures)

    # 删除第一行，即包含'geneSymbol'和'cell_name'的行
    data_frame <- data_frame[-1, ]
    colnames(data_frame) <- c("Gene", "Celltype")

    return(data_frame)
  }

  # 定义将数据框转换为特定格式的函数
  SCINA_markers_transformed <- function(data_frame) {
    data_frame %>%
      group_by(Celltype) %>%
      arrange(Celltype, Gene) %>%
      ungroup()
  }

  # 转换标记基因为数据框
  data_frame <- transform_signature_to_df(signatures)
  result <- SCINA_markers_transformed(data_frame)

  # 输出结果
  print(result)

  # 保存结果到 txt 文件中
  write.table(result, file = "demo.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  marker <- read.table("demo.txt", header = TRUE, sep = "\t")

  ######## Heatmap ######

  #heatmap <- DoHeatmap(scRNA, features = marker$Gene, assay = "RNA")
  #print(heatmap)

  # 取出矩阵并进行标准化
  Avg.Expr <- AverageExpression(scRNA, features = marker[!duplicated(marker$Gene), "Gene"], assay = "RNA")
  Avg.Expr.df <- Avg.Expr$RNA

  # 表达量取log(*+1)处理
  Avg.Expr.df <- log(Avg.Expr.df + 1)

  # 生成 annotation_row
  annotation_row <- marker[!duplicated(marker$Gene), ]
  annotation_row <- as.data.frame(annotation_row$Celltype)
  rownames(annotation_row) <- marker[!duplicated(marker$Gene), "Gene"]

  # 提取列名的数字部分
  col_numbers <- as.numeric(gsub("[^0-9]", "", colnames(Avg.Expr.df)))

  # 获取列的排序索引（从小到大）
  col_order <- order(col_numbers)

  # 使用排序后的索引重新排列表格的列
  Avg.Expr.df <- Avg.Expr.df[, col_order]

  # 生成 annotation_row，但是将名称设置为 "Cell"
  annotation_row <- marker[!duplicated(marker$Gene), ]
  annotation_row <- data.frame(Cell = annotation_row$Celltype, stringsAsFactors = FALSE)
  rownames(annotation_row) <- marker[!duplicated(marker$Gene), "Gene"]

  # 使用pheatmap画热图

  pheatmap(Avg.Expr.df,
           show_colnames = TRUE,
           show_rownames = TRUE,
           border_color = NA,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_row = annotation_row,
           use_raster = FALSE,
           cellwidth = 18,
           cellheight = 6,
           fontsize = 7)

  # 保存热图为PNG文件，使用用户指定的大小和分辨率
  dev.copy(png, file = output_file, width = width, height = height, res = res)
  dev.off()
}
