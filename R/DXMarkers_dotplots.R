#' @title Generate DotPlots from a Seurat object and save as PDF
#' @description This function generates DotPlots for each cell type from a Seurat object, using the genes specified in a csv file, and saves the plots as PDFs in a specified directory.
#' @param pbmc_data a Seurat object
#' @param csv_file a string, path to the csv file containing cell type and gene information
#' @param directory a string, path to the directory where the PDFs will be saved
#' @return NULL
#' @name DXMarkers_dotplots
#' @examples
#' DXMarkers_dotplots(pbmc3k.final, "confirm_SCINA_markers.csv", "C:\\Users\\Wandering\\Desktop\\测试")
#' @export
#' @title Generate DotPlots from a Seurat object and save as PDF
#' @description This function generates DotPlots for each cell type from a Seurat object, using the genes specified in a csv file, and saves the plots as PDFs in a specified directory.
#' @param pbmc_data a Seurat object
#' @param csv_file a string, path to the csv file containing cell type and gene information
#' @param directory a string, path to the directory where the PDFs will be saved
#' @return NULL
#' @examples
#' DXMarkers_dotplots(pbmc3k.final, "confirm_SCINA_markers.csv", "C:\\Users\\Wandering\\Desktop\\测试")
#' @export
library(patchwork)
library(RColorBrewer)
library(grid)

mycolors <- c("#3f5cae","#7bb3d6","#d0ebf3","#fdfac2","#fdc178","#ef6845","#b0182a")

DXMarkers_dotplots <- function(pbmc_data, csv_file, directory) {
  # 从csv文件中读取数据并转换为cell_genes列表
  read_cell_genes <- function(file) {
    data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    cell_genes <- lapply(data, function(x) x[x != ""])
    return(cell_genes)
  }
  cell_genes <- read_cell_genes(csv_file)

  # 创建输出文件夹
  if (!dir.exists(directory)) {
    dir.create(directory)
  }

  # 用于保存所有重复基因及其对应细胞类型的数据框
  duplicated_genes_df <- data.frame(Gene = character(), CellType = character(), stringsAsFactors = FALSE)

  # 批量生成DotPlots
  for (cell_type in names(cell_genes)) {
    genes <- cell_genes[[cell_type]]

    # 检查并保存重复的基因
    dup_genes <- genes[duplicated(genes)]
    if(length(dup_genes) > 0){
      duplicated_genes_df <- rbind(duplicated_genes_df, data.frame(Gene = dup_genes, CellType = cell_type))
    }

    # 移除重复的基因
    genes <- unique(genes)

    p <- DotPlot(pbmc_data, features = genes) +
      coord_flip() +
      ggtitle(cell_type) +
      theme_bw(base_size = 18) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = .5),
            panel.background = element_rect(fill = "white", colour = "black")) +
      labs(x = NULL, y = NULL) +
      guides(size = guide_legend(order = 3)) +
      scale_color_gradientn(values = scales::rescale(seq(0, 1, length.out = length(mycolors))), colours = mycolors) # 修改颜色

    # 保存为PNG文件
    ggsave(file.path(directory, paste0(cell_type, "_DotPlot.png")), p, width = 10, height = 5)
  }

  # 对所有重复基因进行后续分析
  duplicated_genes <- unique(duplicated_genes_df$Gene)
  duplicated_genes <- duplicated_genes[duplicated_genes %in% rownames(pbmc_data)]

  if(length(duplicated_genes) > 0){
    # 创建点图（DotPlot）
    P8 <- DotPlot(pbmc_data, features = duplicated_genes, assay='RNA') + coord_flip() + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black")) +
      scale_color_gradientn(values = scales::rescale(seq(0, 1, length.out = length(mycolors))), colours = mycolors) # 修改颜色

    # 创建密度图（VlnPlot）
    P9 <- VlnPlot(object = pbmc_data, features = duplicated_genes, log = T) + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black"))

    # 创建特征图（FeaturePlot）
    P10 <- FeaturePlot(object = pbmc_data, features = duplicated_genes) + theme_bw() + theme(panel.background = element_rect(fill = "white", colour = "black"))

    # 保存为PNG文件，长宽加大两倍
    ggsave(file.path(directory, "00Duplicated_Genes_DotPlot.png"), P8, width = 8, height = 8)
    ggsave(file.path(directory, "00Duplicated_Genes_VlnPlot.png"), P9, width = 20, height = 20)
    ggsave(file.path(directory, "00Duplicated_Genes_FeaturePlot.png"), P10, width = 20, height = 20)
  }

  # 将重复基因和其对应的细胞类型保存为CSV文件
  write.csv(duplicated_genes_df, file.path(directory, "duplicated_genes.csv"), row.names = FALSE)
}






