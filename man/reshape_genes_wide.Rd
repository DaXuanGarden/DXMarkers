% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape_genes_wide.R
\name{reshape_genes_wide}
\alias{reshape_genes_wide}
\title{重塑基因数据为宽格式并排序，并保存为CSV文件}
\usage{
reshape_genes_wide(sorted_genes, output_dir = "DXMarkers")
}
\arguments{
\item{sorted_genes}{经过排序后的基因数据框}

\item{output_dir}{保存CSV文件的目录。默认为"DXMarkers"文件夹。}
}
\value{
重塑后的基因数据框（宽格式）作为数据框
}
\description{
该函数将经过排序后的基因数据重塑成宽格式，每个cluster为一列，并保存结果至CSV文件中。
}
\examples{
sorted_result <- reshape_genes_wide(top10_data)
}
