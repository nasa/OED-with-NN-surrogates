#' Generate a MAC file
#' @export
genMAC <- function(constituents, analysis_type, loading, failure, name, out_file = paste0(name, ".MAC")){
  library("readr", character.only = TRUE)
  
  constituent_lines <- list("*CONSTITUENTS", NMATS = 2L) 
  list(M=1, CMOD=6, MATID="U", MATDB=1, EL=sprintf("%f,%f,%f,%f,%f,%f,%f"))
  mat_str <- constituent_lines
  readr::write_lines(x = mat_str, path = "test.MAC")
  # NMATS=2
  # # -- Carbon fiber
  # M=1 CMOD=6 MATID=U MATDB=1 &
  #   EL=230995.65703,238326.98225,0.20000,0.20000,97500.00000,0.0E-6,0.0E-6
  # # -- Epoxy matrix
  # M=2 CMOD=6 MATID=U MATDB=1 &
  #   EL=3454.93265,3454.93265,0.35403,0.35403,1275.79385,0.0E-6,0.0E-6
}