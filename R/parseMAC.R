#' parseMac
#' @export
parseMAC <- function(file_path){
  mac_file <- readr::read_lines(file = file_path)
  section_idx <- grep(pattern = "\\*.*", x = mac_file, ignore.case = TRUE)
  con_idx <- grep(pattern = "\\*CONSTITUENTS", x = mac_file, ignore.case = TRUE)
  nmats <- as.integer(gsub(mac_file[con_idx + 1], pattern = "\\s*NMATS=(\\d+)", replacement = "\\1", ignore.case = TRUE))
  ws_num_rgx <- "\\s*([-+]?[0-9]*\\.?[0-9]+(?:[eED][-+]?[0-9]+)?)\\s*"
  con_rgx <- sprintf("M=\\d+\\s*CMOD=\\d+\\s*MATID=\\w\\s*MATDB=\\d+\\s*EL=%s,%s,%s,%s,%s,%s,%s", ws_num_rgx, ws_num_rgx, ws_num_rgx, ws_num_rgx, ws_num_rgx, ws_num_rgx, ws_num_rgx)
  con_props1 <- as.numeric(strsplit(gsub(mac_file[con_idx+2], pattern = con_rgx, replacement = "\\1 \\2 \\3 \\4 \\5 \\6 \\7"), split = "\\s+")[[1]][-1])
  con_props2 <- as.numeric(strsplit(gsub(mac_file[con_idx+3], pattern = con_rgx, replacement = "\\1 \\2 \\3 \\4 \\5 \\6 \\7"), split = "\\s+")[[1]][-1])
  lam_idx <- grep(pattern = "\\*LAMINATE", x = mac_file, ignore.case = TRUE)
  n_ply <- as.numeric(trimws(gsub(x = mac_file[[lam_idx+1L]], pattern = "NLY=(\\d+)", replacement = "\\1", ignore.case = TRUE)))
  lam_data <- mac_file[(lam_idx+2):(lam_idx+n_ply+1)]
  lam_props <- strsplit(lam_data, split = "\\s+")
  lam_data <- as.data.frame(t(sapply(lam_props, function(lp){
    as.numeric(gsub(x = lp[-1], pattern = "(.+)=(.+)", replacement = "\\2"))
  })))
  colnames(lam_data) <- gsub(x = lam_props[[1]][-1], pattern = "(.+)=(.+)", replacement = "\\1")
  
  
  str_lines <- grep(x = mac_file, pattern = sprintf("X11=%s", ws_num_rgx))
  str_props <- as.numeric(gsub(x = mac_file[str_lines], pattern = sprintf(".*X11=(%s).*", ws_num_rgx), replacement = "\\1"))
  con_props1 <- c(con_props1, str_props[[1]])
  con_props2 <- c(con_props2, str_props[[2]])
  
  names(con_props1) <- c("E_a", "E_t", "V_a", "V_t", "G_a", "alpha_a", "alpha_t", "sigma")
  names(con_props2) <- c("E_a", "E_t", "V_a", "V_t", "G_a", "alpha_a", "alpha_t", "sigma")
  
  res <- list(
    fiber_constituent_properties = con_props1, 
    matrix_constituent_properties = con_props2,
    laminate_properties = lam_data
  )
  return(res)
  # mat_props <- getMatches(x = mac_file[[3]], rgx = "M=(\\d+)\\s+CMOD=(\\d+)\\s+MATID=(.*)\\s+MATDB=(\\d+)\\s+EL=(.*)")[[1]]

}


#' writeMAC
#' 
#' @import readr
#' @export
writeMAC <- function(out_path, mat_type, mat1, mat2, lam, mech, sc1, sc2, solver, plot_config, print_level=6L){
  
  ## Constituent definition
  n_mat <- 2L
  constituents_def <- sprintf("*CONSTITUENTS\n  NMATS=%d\n", n_mat)
  mat_def <- " M=%d CMOD=%d MATID=%s MATDB=%d EL=%f,%f,%f,%f,%f,%f,%f\n"
  mat1_def <- paste0("# -- Carbon fiber\n", do.call(sprintf, append(list(mat_def), mat1)))
  mat2_def <- paste0("# -- Epoxy matrix\n", do.call(sprintf, append(list(mat_def), mat2)))
  constituent_conf <- paste0(constituents_def, mat1_def, mat2_def)
  
  ## Laminate definition
  n_ply <- nrow(lam)
  lam_def <- sprintf("*LAMINATE\n  NLY=%d\n", n_ply)
  lam_conf <- paste0(lam_def, parse_df(lam, pre_space = "    "))
  # lam_def <- "LY=%d MOD=%d THK=%f ANG=   %f ARCHID=%d VF=%f F=%d M=%d"
  
  ## Applied Mechanical Loading Definitions
  lop <- 99L # loading option
  mech_def <- sprintf("*MECH\n  LOP=%d\n", lop)
  mech_conf <- paste0(mech_def, parse_df(mech, pre_space = "  "))
    
  ## Subcell static failure section
  sc1_conf <- paste0("  MAT=1 NCRIT=1\n    ", parse_row(sc1), "\n")
  sc2_conf <- paste0("  MAT=2 NCRIT=1\n    ", parse_row(sc2))
  sc_conf <- paste0("*FAILURE_SUBCELL\n", " NMAT=2\n", sc1_conf, sc2_conf)
  
  ## Solver to use 
  solver_conf <- paste0("*SOLVER\n  ", parse_row(solver))
  
  ## Print configuration 
  print_conf <- sprintf("*PRINT\n  NPL=%d", print_level)
  
  ## Output section 
  tags <- strsplit(parse_row(plot_config), split = " ")[[1]]
  tags <- c(append(tags[1:2], paste0(tags[3:5], collapse = " ")), tags[-c(1:5)])
  tags <- paste0(unname(sapply(tags, function(tag) paste0("  ", tag))), collapse = "\n")
  plot_conf <- paste0("*XYPLOT\n", tags)
  
  ## Make the MAC output file
  write_lines(x = c(mat_type, constituent_conf, lam_conf, mech_conf, sc_conf, solver_conf, print_conf, plot_conf, "*END"), path = out_path)  
}

#' parse_df
#' @description Parses each row of a named data.frame into a tag list 
#' @return if collapse=TRUE, a strign of the parsed rows is returned. Otherwise a list where each elements represents a row is returned. 
#' @export
parse_df <- function(df, pre_space = "", collapse=FALSE){
  
  ## Parse row by row
  res <- vector("list", length=nrow(df))
  for (i in seq(nrow(df))){
    res[[i]] <- paste0(pre_space, parse_row(df[i,,drop=FALSE], colnames(df)))
  }
  
  ## Unless collapse is requested, return a list of the parsed rows. Otherwise collapse w/ newline. 
  if (is.logical(collapse) && !collapse){
    return(res)
  } else {
    return(paste0(unlist(res), collapse="\n"))
  }
}

#' parse_row 
#' @description Converts a vector with names into a tag=value space-delimited string
#' @return a character string containing a 'tag=value' for each column of a row.  
#' @export
parse_row <- function(x, x_names = names(x), collapse = " "){
  paste0(as.vector(sapply(x_names, function(h) { 
    paste0(h, "=", switch(readr::guess_parser(x[[h]], guess_integer = TRUE), 
                          "integer"=paste0(as.integer(x[[h]]), collapse = ", "), 
                          "double"=paste0(sprintf("%.5f", x[[h]]), collapse = ", "), 
                          x[[h]]))
  })), collapse = collapse)
}
