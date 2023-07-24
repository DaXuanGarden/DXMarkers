#' @title View Data Source
#' @description This function returns the dataset from a particular data source.
#' @param data_source A character string specifying the data source to view. It can be either "CellMarker", "CellMarker2.0", or "PanglaoDB".
#' @return A dataframe containing the dataset from the selected data source.
#' @examples
#' \dontrun{
#'   view_data <- view_data_source(data_source = "CellMarker")
#' }
view_data_source <- function(data_source) {
  if (data_source == "CellMarker") {
    data <- readRDS(system.file("data", "CellMarker_all_cell_markers.Rds", package = "DXMarkers"))
  } else if (data_source == "CellMarker2.0") {
    data <- readRDS(system.file("data", "CellMarker2_0_Cell_marker_All.Rds", package = "DXMarkers"))
  } else if (data_source == "PanglaoDB") {
    data <- readRDS(system.file("data", "PanglaoDB_markers.Rds", package = "DXMarkers"))
  } else {
    stop("Invalid data source. Please choose from 'CellMarker', 'CellMarker2.0', or 'PanglaoDB'.")
  }
  return(data)
}

