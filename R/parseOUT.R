#' parseOUT
#' 
#' Parses section III of an *.out file produced by MAC/GMC into a matrix 
#' @import readr
#' @import stringr
#' @export
parseOUT <- function(file_path){
  epoxy_file <- read_file(file = file_path)
  out_file <- read_lines(file = file_path)
  idx <- grep(x = out_file, "Exx", ignore.case = FALSE)
  res <- list()
  res[["Exx"]] <- as.numeric(str_sub(out_file[idx], -10, -1))
  res[["Nxy"]] <- as.numeric(str_sub(out_file[idx+1], -10, -1))
  return(res)
}
parseOUT_old <- function(file_path, sections = c(2L, 3L)){
  
  ## Generic regexs
  problem_title_regex <- "\\*+(?:\\s+PROBLEM TITLE\\s+\\*+){5}\\s*([0-9]*(?: \\w+)+)\\s+\\*+\\s+"
  numeric_regex <- "([-+]?[0-9]*\\.?[0-9]+(?:[eED][-+]?[0-9]+)?)"
  ws_numeric_regex <- paste0("\\s*", numeric_regex, "\\s*")
  
  ## Read in the file as a raw character string
  epoxy_file <- read_file(file = file_path)
  out_file <- read_lines(file = file_path)
  sec_idx <- grep(x = out_file, "\\s*\\*+\\s*Section", ignore.case = TRUE)
  section1 <- out_file[sec_idx[1]:(sec_idx[2]-1)]
  section2 <- out_file[sec_idx[2]:(sec_idx[3]-1)]
  section3 <- out_file[sec_idx[3]:length(out_file)]

  ## Simple function to extract the matches
  getMatches <- function(x, rgx){
    rgx_matches <- regmatches(x = x, m = gregexpr(pattern = rgx, text = x, perl = TRUE))
    rgx_tmp <- regmatches(x = rgx_matches[[1]], regexec(pattern = rgx, text = rgx_matches[[1]], perl = TRUE))
    lapply(rgx_tmp, function(m_i) m_i[2:length(m_i)])
  }
  findAndMatch <- function(x, rgx, rmatch = "\\1"){
    gsub(x = grep(pattern = rgx, x = x, value = TRUE), pattern = rgx, replacement = rmatch)
  }
  
  ## The final results
  res <- list()
  
  ## Parse section II 
  if (2L %in% sections){

    toNum <- function(x){ as.numeric(unlist(x)) }
    as_mat <- function(x){ do.call(rbind, lapply(x, as.numeric)) }
    
    ## Get number of layers
    layer_regex <- "\\s*Total Number of Layers\\s*\\(nly\\)\\s*=\\s*(\\d+)"
    n_layers <- as.integer(findAndMatch(section1, rgx = layer_regex))
      
    ## Some properties  properties that  
    layer_vec <- function(x){
      modulo <- length(x)/n_layers
      if (modulo == 1){ return(toNum(x)) } 
      else {
        from <- seq(1, length(x), by = modulo)
        to <- seq(1, length(x), by = modulo) + modulo - 1
        mapply(function(i, j){  toNum(x[i:j]) }, from, to, SIMPLIFY = FALSE)
      }
    }
    
    ## Extract the subcell meta information for each layer 
    layer_start_rgx <- "\\s*\\(BETA\\s*, GAMMA\\)\\s*SUBCELL \\#\\s*SUBCELL MATERIAL\\s*SUBCELL VOLUME\\s*"
    layer_end_rgx <- paste0("\\s*TOTAL VOLUME =\\s*", numeric_regex)
    idx_start <- gregexpr(pattern = layer_start_rgx, text = epoxy_file, perl = TRUE)
    idx_end <- gregexpr(pattern = layer_end_rgx, text = epoxy_file, perl = TRUE)
    subcell_info_idx <- cbind(as.vector(idx_start[[1]])+attr(idx_start[[1]], "match.length"), 
                              as.vector(idx_end[[1]]))
    subcell_meta_rgx <- paste0("(\\d+)\\s*,\\s*(\\d+)\\s*(\\d+)\\s*(\\d+)\\s*", numeric_regex)
    subcell_meta <- lapply(1:nrow(subcell_info_idx), function(i){
      idx <- subcell_info_idx[i,]
      subcell_data <- read_lines(substr(epoxy_file, start = idx[1], stop = idx[2]))
      tmp <- lapply(subcell_data, function(sm) as.numeric(getMatches(sm, subcell_meta_rgx)[[1]]))
      res <- do.call(rbind, tmp)
      colnames(res) <- c("beta", "gamma", "subcell_num", "subcell_material", "subcell_volume")
      return(res)
    })
    subcell_meta <- data.table::rbindlist(lapply(subcell_meta, as.data.frame), idcol = "layer")
    
    ## Extract the scalar statistics associated with each layer 
    total_volume <- as_mat(getMatches(epoxy_file, layer_end_rgx))
    material_volume_ratio <- as_mat(layer_vec(getMatches(epoxy_file, paste0("MATERIAL NO.=\\s*\\d+\\s+VOLUME RATIO=\\s*", numeric_regex))))
    RUC_dim <- as_mat(getMatches(epoxy_file, sprintf("RUC Dimensions: H = %s\\s*L = %s", numeric_regex, numeric_regex)))
    subcell_dim <- as_mat(getMatches(epoxy_file, sprintf("Subcell Dimensions:\\s*h =\\s*%s %s\\s*l =\\s*%s %s\\s*", numeric_regex, numeric_regex, numeric_regex, numeric_regex)))
    
    ## Get the CI/CG Effective/Macro Stiff Matrix information 
    CG_mat <- as_mat(getMatches(x = epoxy_file, paste0("CG - Effective/Macro Stiffness Matrix\\s*", paste0(rep(ws_numeric_regex, 12L), collapse = ""))))
    CI_mat <- as_mat(getMatches(x = epoxy_file, paste0("CI - Effective/Macro Compliance Matrix\\s*", paste0(rep(ws_numeric_regex, 12L), collapse = ""))))
    
    ## Get the engineering moduli
    engineering_moduli_regex <- sprintf("Effective Engineering Moduli\\s*E11S=%sN12S=%sN13S=%sE22S=%sN23S=%sE33S=%sG23S=%sG13S=%sG12S=%s", ws_numeric_regex, ws_numeric_regex, ws_numeric_regex, ws_numeric_regex, ws_numeric_regex, ws_numeric_regex, ws_numeric_regex, ws_numeric_regex, ws_numeric_regex) 
    eng_moduli <- as_mat(getMatches(epoxy_file, engineering_moduli_regex))
    
    ## Get the thermal expansion coefficients
    exp_coeff <- as_mat(getMatches(epoxy_file, sprintf("Effective Tangent Thermal Expansion Coefficients\\s*%s%s%s", ws_numeric_regex, ws_numeric_regex, ws_numeric_regex)))
    
    ## Local Q Stiffness layer
    local_q_stiffness_rgx <- paste0("Local Q Stiffness For Layer\\s*\\d+\\s*Integration Pt.\\s*\\d+\\s*", paste0(rep(ws_numeric_regex, 9L), collapse = ""))
    local_q_stiffness <- as_mat(layer_vec(getMatches(epoxy_file, rgx = local_q_stiffness_rgx)))
    
    ## Global Q Stiffness layer
    global_q_stiffness_rgx <- paste0("Global Q Stiffness For Layer\\s*\\d+\\s*Integration Pt.\\s*\\d+\\s*", paste0(rep(ws_numeric_regex, 9L), collapse = ""))
    global_q_stiffness <- as_mat(layer_vec(getMatches(epoxy_file, rgx = global_q_stiffness_rgx)))
    
    ## Save the results
    tmp_res <- cbind(total_volume, material_volume_ratio, RUC_dim, subcell_dim, CG_mat, CI_mat, 
                     eng_moduli, exp_coeff, local_q_stiffness, global_q_stiffness)
    colnames(tmp_res) <- c("total_volume", paste0("volume_ratio_", 1:ncol(material_volume_ratio)), 
                           "ruc_h", "ruc_l", paste0("subcell_dim_", c("h1", "l1", "h2", "l2")), 
                           paste0("CG_", 1:ncol(CG_mat)), paste0("CI_", 1:ncol(CI_mat)), 
                           "E11S", "N12S", "N13S", "E22S", "N23S", "E33S", "G23S", "G13S", "G12S", 
                           paste0("tan_therm_exp_", 1:3), 
                           paste0("local_q_stiffness_", 1:ncol(local_q_stiffness)),
                           paste0("global_q_stiffness_", 1:ncol(global_q_stiffness))) 
    res[["section II"]][["layer_statistics"]] <- tmp_res
    
    ## Get initial ABD matrix 
    a_idx <- grep("\\s*Laminate Axial Stiffness Matrix \\[A\\]", x = section2, ignore.case = TRUE)
    b_idx <- grep("\\s*Laminate Coupling Stiffness Matrix \\[B\\]", x = section2, ignore.case = TRUE)
    d_idx <- grep("\\s*Laminate Bending Stiffness Matrix \\[D\\]", x = section2, ignore.case = TRUE)
    A <- as.numeric(strsplit(paste0(section2[(a_idx+2):(a_idx+4)], collapse = ""), split = "\\s+")[[1]][-1])
    B <- as.numeric(strsplit(paste0(section2[(b_idx+2):(b_idx+4)], collapse = ""), split = "\\s+")[[1]][-1])
    D <- as.numeric(strsplit(paste0(section2[(d_idx+2):(d_idx+4)], collapse = ""), split = "\\s+")[[1]][-1])
    res[["section II"]][["ABD_matrix"]] <- list(A, B, D)
  }
  
  ## Parse section III 
  if (3 %in% sections){
    ## Regex's that will be needed
    time_temp_tstep_regex <- paste0("\\s*\\d\\s*TIME:\\s*", numeric_regex, "\\s*TEMP:\\s*", numeric_regex, "\\s*TSTEP:\\s*", numeric_regex, "\\s*")
    force_moment_regex  <- paste0("\\s*FORCE\\(N\\), MOMENT\\(M\\):", paste0(rep(ws_numeric_regex, 6L), collapse = ""))
    strain_curvature_regex  <- paste0("\\s*STRAIN, CURVATURE:", paste0(rep(ws_numeric_regex, 6L), collapse = ""))
    inelastic_nm_regex  <- paste0("\\s*INELASTIC N, M:", paste0(rep(ws_numeric_regex, 6L), collapse = ""))
    out_of_plane_strain_regex  <- paste0("\\s*OUT-OF-PLANE STRAIN:", paste0(rep(ws_numeric_regex, 3L), collapse = ""))
    adb_row_regex <- paste0("\\s*\\|", paste0(rep(ws_numeric_regex, 6L), collapse = ""), "\\s*\\|\\s*", collapse = "")
    abd_matrix_regex  <- paste0("ABD MATRIX:", paste0(rep(adb_row_regex, 6L), collapse = ""), collapse = "")
    
    
    ## Given a source text 'x' and a regex 'rgx', extracts an (nx[n_numeric]) matrix of detected lines of numbers
    extractNumericRegex <- function(x, n_numeric, rgx){
      rgx_matches <- regmatches(x = x, m = gregexpr(pattern = rgx, text = x, perl = TRUE))
      rgx_tmp <- regmatches(x = rgx_matches[[1]], regexec(pattern = rgx, text = rgx_matches[[1]], perl = TRUE))
      do.call(rbind, lapply(rgx_tmp, function(x) as.numeric(stringr::str_replace(x[2:(n_numeric+1)], "D", "e"))))
    }
    
    ## Extract the: 
    ## 1. time, temperature, and T-step
    ## 2. force and moment vectors
    ## 3. strain and curvature vectors
    ## 4. inelastic N,M vectors 
    ## 5. out-of-plane strain vectors
    time_temp_tstep <- extractNumericRegex(x = epoxy_file, n_numeric = 3, rgx = time_temp_tstep_regex)
    fm <- extractNumericRegex(x = epoxy_file, n_numeric = 6, rgx = force_moment_regex)
    strain_curvature <- extractNumericRegex(x = epoxy_file, n_numeric = 6, rgx = strain_curvature_regex)
    inelastic_nm <- extractNumericRegex(x = epoxy_file, n_numeric = 6, rgx = inelastic_nm_regex)
    out_of_plane_strain <- extractNumericRegex(x = epoxy_file, n_numeric = 3, rgx = out_of_plane_strain_regex)
    
    ## Extract the ABD matrix 
    abd_matches <- regmatches(x = epoxy_file, m = gregexpr(pattern = abd_matrix_regex, text = epoxy_file, perl = TRUE))
    abd_tmp <- regmatches(x = abd_matches[[1]], regexec(pattern = abd_matrix_regex, text = abd_matches[[1]], perl = TRUE))
    abd_matrices <- lapply(abd_tmp, function(x) matrix(as.numeric(stringr::str_replace(x[2:37], "D", "e")), nrow = 6, ncol = 6))
    
    ## Put it all together
    laminate_info <- cbind(time_temp_tstep, fm, strain_curvature, inelastic_nm, out_of_plane_strain)
    XYZ <- function(str){ paste(str, c("X", "Y", "Z"), sep = "")}
    colnames(laminate_info) <- c("time", "temp", "tstep", 
                                 XYZ("force"), XYZ("moment"),
                                 XYZ("strain"), XYZ("curvature"), 
                                 XYZ("inelasticN_"), XYZ("inelasticM_"), 
                                 XYZ("out_of_plane_strain"))
    res[["section III"]] <- laminate_info
    res[["ABD_matrices"]] <- abd_matrices
  }
  
  
  ## Return the results
  return(res)
}