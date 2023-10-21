#' DXChat_group: a function for analyzing and saving data for multiple groups in scRNA data
#'
#' This function takes a Seurat object with RNA counts and metadata, a directory to save output,
#' and a CellChatDB, and performs CellChat analysis for each group specified in the metadata.
#' Results are saved to the output directory with filenames containing the group names.
#'
#' @param scRNA A Seurat object containing scRNA-seq data
#' @param output_dir A string specifying the directory to save output files
#' @param CellChatDB The CellChat database to use for the analysis
#' @return NULL
#' @export
#' @examples
#' # Assuming you have a Seurat object 'seurat' and a CellChat database 'CellChatDB'
#' DXChat_group(seurat, "path/to/output/directory", CellChatDB)
DXChat_group <- function(scRNA, output_dir, CellChatDB) {
  library(CellChat)
  library(patchwork)
  library(ggalluvial)
  library(igraph)
  library(dplyr)

  # 获取所有独特的分组
  groups <- unique(scRNA@meta.data$disease)

  # 为每个独特的分组运行分析
  for (group in groups) {
    print(paste0("Analyzing group ", group, "..."))

    # 提取数据和元数据信息
    data.input <- scRNA[["RNA"]]@counts
    meta <- scRNA@meta.data

    # 提取满足条件的meta信息
    cell.use <- rownames(meta)[meta$disease == group]
    data.input <- data.input[, cell.use]
    meta <- meta[cell.use, ]

    # 删除无效因子
    meta$celltype <- droplevels(meta$celltype)

    # 创建cellchat对象，用meta中的注释作为分组依据
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "celltype")

    # 将选定的数据库内容载入cellchat对象
    cellchat@DB <- CellChatDB

    # 数据预处理
    cellchat <- subsetData(cellchat)
    cellchat@idents <- droplevels(cellchat@idents)

    # 寻找高表达的基因和通讯路径
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # 投影到PPI (Protein-Protein Interaction)
    cellchat <- projectData(cellchat, PPI.human)

    # 计算通讯概率
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = F)

    # 提取细胞间的通讯预测结果
    df.net <- subsetCommunication(cellchat)

    # 计算每对配受体之间的通讯概率；结果保存在"net"插槽中
    cellchat <- computeCommunProbPathway(cellchat)

    # 结果在"netP"插槽中保存了每条通讯路径的预测结果
    cellchat <- aggregateNet(cellchat)

    # 保存结果
    saveRDS(cellchat, file.path(output_dir, paste0(group, ".cellchat.rds")))
    write.csv(df.net, file.path(output_dir, paste0(group, ".df.net.csv")))
  }

  print("All groups analyzed.")
}


