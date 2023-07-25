#' @title 重塑基因数据为宽格式并排序，并保存为CSV文件
#' @description 该函数将经过排序后的基因数据重塑成宽格式，每个cluster为一列，并保存结果至CSV文件中。
#' @param sorted_genes 经过排序后的基因数据框
#' @return 一个包含重塑后的基因数据框的list
#' @export
#' @examples
reshape_genes_wide <- function(sorted_genes) {
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

  # 创建DXMarkers文件夹（如果不存在）
  if (!dir.exists("DXMarkers")) {
    dir.create("DXMarkers")
  }

  # 输出排序后的结果并保存为CSV文件
  write.csv(genes_wide, file = "DXMarkers/sorted_genes_wide.csv", row.names = FALSE)

  return(list(sorted_genes_wide = genes_wide))
}
