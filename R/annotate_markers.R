#' @title Annotate Markers
#' @description This function filters a preloaded dataset based on provided top 10 genes for each cluster, species, and tissue class. It also allows to view the data from a particular data source.
#' @param top10_genes_data A named list where names are cluster identifiers and values are vectors of top 10 gene symbols for each cluster.
#' @param species A character string specifying the species. It can be either "Human" or "Mouse".
#' @param tissue_class A character string specifying the tissue class.
#' @param data_source A character string specifying the data source to view. It can be either "CellMarker", "CellMarker2.0", or "PanglaoDB".
#' @return A list with two components: "all_filtered_data" which is a list where each element is a data frame of filtered data for each cluster, and "all_filtered_data_df" which is a data frame combining all the filtered data.
#' @export
#' @examples
#' \dontrun{
#'   data <- list("cluster1" = c("gene1", "gene2"), "cluster2" = c("gene3", "gene4"))
#'   result <- annotate_markers(data, "Human", "Brain")
#'   view_data <- annotate_markers(data_source = "CellMarker")
#' }
annotate_markers <- function(top10_genes_data = NULL, species = NULL, tissue_class = NULL, data_source = NULL) {
  # Load data
  DXMarkers_data <- readRDS(system.file("data", "DXMarkers_data.Rds", package = "DXMarkers"))

  if (!is.null(data_source)) {
    # Subset the dataset based on the data source
    DXMarkers_data <- DXMarkers_data %>% filter(source == data_source)

    # Create the DXMarkers folder if it doesn't exist
    if (!dir.exists("DXMarkers")) {
      dir.create("DXMarkers")
    }

    # Save the filtered data as CSV in the DXMarkers folder
    write.csv(DXMarkers_data, file = "DXMarkers/filtered_data.csv", row.names = FALSE)

    return(DXMarkers_data)
  }

  # Define mapping of species identifiers between data sets
  species_map <- list("Mouse" = c("Mouse", "Mm"), "Human" = c("Human", "Hs"))

  # Create a dataframe with the appropriate species and tissue class
  species_tissue_df <- data.frame(species = species_map[[species]], tissue_class = tissue_class)

  # Subset the dataset based on species and tissue class
  DXMarkers_data <- DXMarkers_data %>%
    semi_join(species_tissue_df, by = c("species", "tissue_class"))

  # Ensure all PMID are of the same type (character)
  DXMarkers_data$PMID <- as.character(DXMarkers_data$PMID)

  # Initialize an empty list to store results
  all_filtered_data <- list()

  # Loop over the clusters in the input data
  for (cluster in names(top10_genes_data)) {
    # Define the gene symbol list for this cluster
    gene_list <- top10_genes_data[[cluster]]

    # Filter the data by the gene symbol list
    filtered_data <- DXMarkers_data %>% filter(geneSymbol %in% gene_list)

    # Only add the cluster as a new column if the filtered data is not empty
    if (nrow(filtered_data) > 0) {
      # Add the cluster as a new column in the filtered data
      filtered_data$cluster <- cluster

      # Add the filtered data to the list of results
      all_filtered_data[[cluster]] <- filtered_data
    }
  }

  # Combine all the filtered data into one data frame
  all_filtered_data_df <- bind_rows(all_filtered_data, .id = "cluster")

  # Create the DXMarkers folder if it doesn't exist
  if (!dir.exists("DXMarkers")) {
    dir.create("DXMarkers")
  }

  # Save the combined data frame as CSV in the DXMarkers folder
  write.csv(all_filtered_data_df, file = "DXMarkers/combined_filtered_data.csv", row.names = FALSE)

  # Return the combined data frame and the list of results
  return(list(all_filtered_data = all_filtered_data, all_filtered_data_df = all_filtered_data_df))
}
