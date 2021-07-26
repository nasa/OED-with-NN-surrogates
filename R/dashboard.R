#' dashboard
#' @description Runs the Material Analysis dashboard.
#' @param X material dataset (optional)
#' @import data.table shiny
#' @export
dashboard <- function(X = NULL, dev_mode = TRUE){
  if (!missing(X)){
    dimensions <- colnames(X)
    X <- parseData(X, dim = dimensions, all_numeric = FALSE, remove.NA = FALSE, scale = FALSE)
    X <- X[, lapply(.SD, function(x_i) if (is.character(x_i) || is.factor(x_i)) as.factor(x_i) else as.numeric(x_i))]
  }
  
  ## If developing the app, set to true to work from the files directly. FALSE will only use the files loaded with 
  ## the MaterialsAnalysis package. 
  if (dev_mode){
    app_dir <- file.path("/home", "mpiekenb", "grc-materials-analysis", "MaterialAnalysis", "inst", "material_dashboard")
    message(sprintf("Development mode\nRunning app with files at: %s", app_dir))
  } else {
    app_dir <- system.file(package = "MaterialAnalysis")
  }
  shiny::runApp(appDir = app_dir)
  
  # dashboard_file <- system.file(file.path('material_dashboard', 'ui_components', 'psp_dashboard_ui.R'), package = "MaterialAnalysis") 
  # server_file <- system.file(file.path('material_dashboard', 'reactive', 'psp_dashboard_server.R'), package = "MaterialAnalysis") 
  # source(dashboard_file, local = TRUE)
  # source(server_file, local = TRUE)
  # shiny::runApp(shinyApp( ui = ui, server = server ))
}