#' @title Extract Marker Genes for Each Cell Type
#' @name SCINA_markers
#' @description This function extracts the marker genes for each cell type from the result of the 'annotate_markers' function.
#' @param annotate_markers_result The result of the 'annotate_markers' function.
#' @return A data frame where each column represents a cell type and each row in that column represents a gene. The name of the column is the cell type and the entries in that column are the genes.
#' @export
#' @examples
#' \dontrun{
#'   data <- list("cluster1" = c("gene1", "gene2"), "cluster2" = c("gene3", "gene4"))
#'   annotate_result <- annotate_markers(data, "Human", "Brain")
#'   markers <- SCINA_markers(annotate_result)
#' }
library(dplyr)
library(tidyr)
SCINA_markers <- function(annotate_markers_result, output_file = "SCINA_markers.csv") {

  cell_type_genes_df <- annotate_markers_result %>%
    select(cell_name, geneSymbol) %>%   # 提取cell_name 和 geneSymbol
    group_by(cell_name) %>%             # 按照cell_name 进行分组
    mutate(id = row_number()) %>%        # 创建一个辅助列，用于后续的分布
    pivot_wider(names_from = cell_name, values_from = geneSymbol) %>% # 展开表格
    select(-id) %>%                     # 删除id列
    replace(is.na(.), "")               # 替换NA值为 ""

  write.csv(cell_type_genes_df, file = output_file, row.names = FALSE, na = "")

  return(cell_type_genes_df)
}

