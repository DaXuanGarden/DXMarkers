#' Seurat_to_matrix
#'
#' The function is designed to extract the gene expression matrix from a Seurat object.
#' It will create a new matrix where each row represents a gene, and each column represents a cell.
#' Note: The input should be a Seurat object.
#'
#' @param seurat_obj A Seurat object.
#'
#' @return A gene expression matrix where each row represents a gene and each column represents a cell.
#'
#' @examples
#' \dontrun{
#' seurat_data <- Read10X(data.dir = "../data")
#' seurat_obj <- CreateSeuratObject(counts = seurat_data, project = "seurat_data", min.cells = 3, min.features = 200)
#' expression_matrix <- Seurat_to_matrix(seurat_obj)
#' }
#' @export
Seurat_to_matrix <- function(seurat_obj) {
  # Extract the count matrix and convert to a data frame
  expression_matrix <- as.data.frame(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))

  # Write to a temporary file
  temp_file <- tempfile(fileext = ".csv")
  data.table::fwrite(x = expression_matrix, row.names = T, file = temp_file)

  # Read back into a matrix
  expression_matrix <- as.matrix(read.csv(temp_file, row.names = 1))

  return(expression_matrix)
}
