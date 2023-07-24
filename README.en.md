# DXMarkers: A One-Click Solution for Retrieving Marker Genes from "CellMarker", "CellMarker2.0", and "PanglaoDB" Databases

In single-cell RNA sequencing (scRNA-seq) analysis, distinguishing cell types is a crucial step. This process often relies on marker genes, which are highly expressed in certain cell types. DXMarkers, a well-designed R package, is designed to achieve this goal. It can efficiently retrieve specific genes from three major marker gene databases: "CellMarker", "CellMarker2.0", and "PanglaoDB", and list the corresponding cell types. Additionally, DXMarkers supports filtering search results by species and tissue type to provide more precise information.

## Installation Guide

First, you need to install the DXMarkers package in your R environment. Below are the steps to install the package locally:

```R
# Set the working path
setwd("/home/data/t050446/01 Single Cell Project/Acute Pancreatitis")

# Load the devtools package for local R package installation
library(devtools)

# Install the local version of the DXMarkers package
install_local("/home/data/t050446/Non-tumor research/DXMarkers_1.0.tar.gz")
```


## One-Click Retrieval of Marker Genes

After performing scRNA-seq data analysis and identifying the marker genes of each cell population, you can use the `annotate_markers` function of DXMarkers to search the above three databases at once and get the cell type corresponding to each gene:

```R
# Load the DXMarkers package
library(DXMarkers)

# Read the marker gene data
top10_genes_data <- read.csv("top10_sorted_genes.csv", stringsAsFactors = FALSE)

# Annotate genes using the annotate_markers function, using the mouse pancreas as an example
DXMarkers_result <- annotate_markers(top10_genes_data, "Mouse", "Pancreas")
```


## Access Built-In Data Sources

DXMarkers also provides a feature to view built-in data sources, allowing you to browse the complete datasets in "CellMarker", "CellMarker2.0", or "PanglaoDB":

```R
# Load the DXMarkers package
library(DXMarkers)

# Specify the data source
data_source <- "CellMarker"

# Browse the data of the specified data source
data <- view_data_source(data_source = data_source)
```


DXMarkers is designed to assist users in finding key marker genes in large-scale single-cell RNA-seq data and quickly and accurately retrieve corresponding information from various databases. Whether you are a beginner in single-cell data analysis or an experienced researcher, DXMarkers can help you save time and improve work efficiency.

