DXMarkers_check <- function(scRNA, output_png = "heatmap_to_check.png", width = 7, height = 5, dpi = 300) {
  # 载入所需的包
  library(Seurat)
  library(RColorBrewer)
  library(qs)
  library(ComplexHeatmap)
  library(dplyr)

  # 创建颜色映射并保存到文件中
  color_map <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(scRNA@active.ident)))
  names(color_map) <- levels(scRNA@active.ident)
  qsave(color_map, 'jdb_color_maps.qs')
  jdb_color_maps <- qread('jdb_color_maps.qs')

  # 找到所有标记基因并保存到文件中
  markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  qsave(markers, 'markers_res_1.5.qs')

  # 读取 markers
  markers <- qread('markers_res_1.5.qs')

  # 选取每个 cluster 中 avg_log2FC 最大的前 10 个基因
  first_marker_list <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

  DoHeatmapPlot <- function(object, groupBy, features) {
    # 载入所需的包
    require(ComplexHeatmap)

    # 获取数据
    plot_data <- SeuratObject::FetchData(object = object,
                                         vars = c(features, groupBy),
                                         slot = 'counts') %>%
      dplyr::mutate(across(.cols = where(is.numeric), .fns = ~ log(.x + 1))) %>%
      dplyr::rename(group = groupBy) %>%
      dplyr::arrange(group)

    clusterInfo <- plot_data$group
    plot_data <- plot_data %>% dplyr::select(-group) %>% t()

    # 创建颜色映射
    col <- jdb_color_maps[1:length(unique(clusterInfo))]
    names(col) <- as.character(unique(clusterInfo))

    # 创建列（细胞类型）注释
    top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col),
                                                       labels = as.character(unique(clusterInfo)),
                                                       labels_gp = gpar(cex = 1.5,
                                                                        col = 'white',
                                                                        family = 'Arial',
                                                                        rot = 45)))

    # 绘制热图
    ht <- Heatmap(matrix = plot_data,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  show_column_names = FALSE,
                  show_row_names = TRUE, #保留基因名
                  show_heatmap_legend = TRUE,
                  column_split = clusterInfo,
                  top_annotation = top_anno,
                  column_title = NULL,
                  use_raster = FALSE,
                  heatmap_legend_param = list(
                    title = 'log(count+1)',
                    title_position = 'leftcenter-rot'
                  ))

    # 返回 ht 对象
    return(ht)
  }
  # 生成 PNG 文件
  png(file = output_png, width = width, height = height, units = "in", res = dpi)
  heatmap_plot <- DoHeatmapPlot(object = scRNA, groupBy = 'celltype', features = unique(unlist(unique(first_marker_list))))
  print(heatmap_plot)
  dev.off()

  # 返回 PNG 图像的文件路径
  return(output_png)
}

