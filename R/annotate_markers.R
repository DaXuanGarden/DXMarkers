#' Annotate Cell and Gene Markers from Multiple Databases
#'
#' This function annotates cell and gene markers from multiple databases.
#'
#' @param top10genes_path file path of the top 10 genes csv file
#' @param species species name, 'Mouse' or 'Human'
#' @param tissue_class tissue class to filter
#'
#' @return a dataframe of the filtered data
#' @importFrom dplyr filter select mutate case_when semi_join
#' @importFrom tidyr separate_rows
#' @importFrom readr read_csv
#' @importFrom magrittr %>%
#'
#' @export
annotate_markers <- function(top10_genes_data, species, tissue_class) {
  # Load data
  DXMarkers_data <- readRDS(system.file("data", "DXMarkers_data.Rds", package = "DXMarkers"))

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
  for(cluster in names(top10_genes_data)) {
    # Define the gene symbol list for this cluster
    gene_list <- top10_genes_data[[cluster]]

    # Filter the data by the gene symbol list
    filtered_data <- DXMarkers_data %>% filter(geneSymbol %in% gene_list)

    # Only add the cluster as a new column if the filtered data is not empty
    if(nrow(filtered_data) > 0) {
      # Add the cluster as a new column in the filtered data
      filtered_data$cluster <- cluster

      # Add the filtered data to the list of results
      all_filtered_data[[cluster]] <- filtered_data
    }
  }

  # Combine all the filtered data into one data frame
  all_filtered_data_df <- bind_rows(all_filtered_data, .id = "cluster")

  # Rename the 'source' column to 'DataSource'
  all_filtered_data_df <- rename(all_filtered_data_df, DataSource = source)

  # Return the combined data frame and the list of results
  return(list(all_filtered_data = all_filtered_data, all_filtered_data_df = all_filtered_data_df))
}
