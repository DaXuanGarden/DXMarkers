#' @title 重塑基因数据为宽格式并排序，并保存为CSV文件
#' @description 该函数将经过排序后的基因数据重塑成宽格式，每个cluster为一列，并保存结果至CSV文件中。
#' @param sorted_genes 经过排序后的基因数据框
#' @param output_dir 保存CSV文件的目录。默认为"DXMarkers"文件夹。
#' @return 重塑后的基因数据框（宽格式）作为数据框
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr arrange group_by mutate
#' @examples
#' sorted_result <- reshape_genes_wide(top10_data)
#' @export
reshape_genes_wide <- function(sorted_genes, output_dir = "DXMarkers") {
  # 检查输入数据框的列是否正确
  required_cols <- c("cluster", "avg_log2FC", "gene")
  if (!all(required_cols %in% colnames(sorted_genes))) {
    stop("输入数据框必须包含以下列: ", paste0("'", paste(required_cols, collapse = "', '"), "'"))
  }

  # 将数据按avg_log2FC从大到小排序，然后重塑成宽格式，每个cluster为一列
  genes_wide <- sorted_genes %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    group_by(cluster) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = cluster, values_from = gene, id_cols = row, values_fill = "")

  # 删除辅助的行编号列
  genes_wide$row <- NULL

  # 修改列名，给每个cluster列名前加上"cluster"
  colnames(genes_wide)[-1] <- paste0("cluster", colnames(genes_wide)[-1])

  # 修改第一列名为"cluster0"
  colnames(genes_wide)[1] <- "cluster0"

  # 创建输出目录（如果不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # 输出排序后的结果并保存为CSV文件
  output_file <- file.path(output_dir, "sorted_genes_wide.csv")
  write.csv(genes_wide, file = output_file, row.names = FALSE)

  return(genes_wide)
}
