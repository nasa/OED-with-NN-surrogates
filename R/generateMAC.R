#' Generates a MAC file
#' @description This function generates a MAC file for simulating tensile tests with MAC/GMC. 
#' @import purrr
#' @param constituents valid constituent data.frame. 
#' @param analysis_type valid constituent data.frame. 
#' @param loading valid constituent data.frame. 
#' @param failure valid constituent data.frame. 
#' @export
generateMAC <- function(folder,
                        number_to_generate_per_laminate_type = 500,
                        delete_old = TRUE,
                        uniform_variation = TRUE,
                        n_ply = 24,
                        # Fiber properties
                        E_fiber_nominal = 338200,
                        E_fiber_percent_variation = 3,
                        UTS_fiber_nominal = 3500,
                        UTS_fiber_percent_variation = 5,
                        nu_fiber_nominal = 0.3,
                        nu_fiber_percent_variation = 3,
                        # Matrix properties
                        E_matrix_nominal = 3450,
                        E_matrix_percent_variation = 3,
                        UTS_matrix_nominal = 80,
                        UTS_matrix_percent_variation = 5,
                        nu_matrix_nominal = 0.35,
                        nu_matrix_percent_variation = 3,
                        # Laminate properties
                        volume_fraction_nominal = 0.6,
                        volume_fraction_variation = 0.03,
                        angle_degree_variation = 1,
                        max_strain=0.2){
  
  # Create directory if it doesn't exist
  if (!file.exists(folder)){
    message(sprintf("Creating new folder: %s", folder))
    dir.create(folder, recursive=TRUE)
  }
  
  # Delete old files
  if (delete_old) {
    message("Deleting old files in input folder...")
    sapply(list.files(folder), function(fn) file.remove(file.path(folder, fn)))
  }
  
  props <- list(
    E_fiber_nominal = E_fiber_nominal,
    UTS_fiber_nominal = UTS_fiber_nominal,
    nu_fiber_nominal = nu_fiber_nominal,
    E_matrix_nominal = E_matrix_nominal,
    UTS_matrix_nominal = UTS_matrix_nominal,
    nu_matrix_nominal = nu_matrix_nominal,
    volume_fraction_nominal = volume_fraction_nominal
  )
  
  # print(E_fiber_nominal)
  # print(UTS_fiber_nominal)
  # print(nu_fiber_nominal)
  # print(E_matrix_nominal)
  # print(UTS_matrix_nominal)
  # print(nu_matrix_nominal)
  # print(volume_fraction_nominal)
  
  props = cross(props)
  
  for (i in 1:length(props)){
    props[[i]]$number_to_generate_per_laminate_type = number_to_generate_per_laminate_type
    props[[i]]$folder = folder
    props[[i]]$start_index = (i - 1) * number_to_generate_per_laminate_type*length(layups) + 1
    props[[i]]$n_ply = n_ply
    props[[i]]$E_fiber_percent_variation = E_fiber_percent_variation
    props[[i]]$UTS_fiber_percent_variation = UTS_fiber_percent_variation
    props[[i]]$nu_fiber_percent_variation = nu_fiber_percent_variation
    props[[i]]$E_matrix_percent_variation = E_matrix_percent_variation
    props[[i]]$nu_matrix_percent_variation = nu_matrix_percent_variation
    props[[i]]$UTS_matrix_percent_variation = UTS_matrix_percent_variation
    props[[i]]$volume_fraction_variation = volume_fraction_variation
    props[[i]]$angle_degree_variation = angle_degree_variation
    props[[i]]$uniform_variation = uniform_variation
    props[[i]]$max_strain = max_strain
  }
  map(props, lift(genOneSetMAC))
  
}

genOneSetMAC <- function(folder, 
              start_index,
              number_to_generate_per_laminate_type,
              uniform_variation,
              n_ply,
              # Fiber properties
              E_fiber_nominal,
              E_fiber_percent_variation,
              UTS_fiber_nominal,
              UTS_fiber_percent_variation,
              nu_fiber_nominal,
              nu_fiber_percent_variation,
              # Matrix properties
              E_matrix_nominal,
              E_matrix_percent_variation,
              UTS_matrix_nominal,
              UTS_matrix_percent_variation,
              nu_matrix_nominal,
              nu_matrix_percent_variation,
              # Laminate properties
              volume_fraction_nominal,
              volume_fraction_variation,
              angle_degree_variation,
              max_strain=0.2){
  
  truncate_between <- function(x, a = 0, b = 1){
    x[x < a] <- a
    x[x > b] <- b
    return(x)
  }
  
  for (li in 1:length(layups)){ assign(paste0("lam", li), value = makeLayup(layups[[li]])) }
  lam_layups <- sapply(1:length(layups), function(i) eval(as.symbol(sprintf("lam%d", i))))
  
  # Simulating graphite and epxoy laminates
  cc <- start_index
  for (i in 1:number_to_generate_per_laminate_type){
    
    # Calculate fiber realization
    E_fiber = property_realization(E_fiber_nominal, E_fiber_percent_variation, uniform=uniform_variation)
    UTS_fiber = property_realization(UTS_fiber_nominal, UTS_fiber_percent_variation, uniform=uniform_variation)
    nu_fiber = property_realization(nu_fiber_nominal, nu_fiber_percent_variation, uniform=uniform_variation)
    #G_fiber = E_fiber / (2 + (1 + nu_fiber))
    G_fiber = E_fiber / (2 * (1 + nu_fiber))
    mat1 <- list(M=1, CMOD=6, MATID="U", MATDB=1)
    mat1 <- append(mat1, list(EL=c(E_fiber,E_fiber,nu_fiber,nu_fiber,G_fiber,-0.68,9.74))) # moduli in MPa
    
    
    # Calculate matrix realization
    E_matrix = property_realization(E_matrix_nominal, E_matrix_percent_variation, uniform=uniform_variation)
    #UTS_matrix = property_realization(UTS_matrix_nominal, UTS_matrix_percent_variation, uniform=uniform_variation)
    nu_matrix = property_realization(nu_matrix_nominal, nu_matrix_percent_variation, uniform=uniform_variation)
    #G_matrix = E_matrix / (2 + (1 + nu_matrix))
    G_matrix = E_matrix / (2 * (1 + nu_matrix))
    UTS_matrix = property_realization(UTS_matrix_nominal, UTS_matrix_percent_variation, uniform=uniform_variation)
    mat2 <- list(M=2, CMOD=6, MATID="U", MATDB=1)
    mat2 <- append(mat2, list(EL=c(E_matrix,E_matrix,nu_matrix,nu_matrix,G_matrix,45,45))) # moduli in MPa
    UTS_matrix_shear <- UTS_matrix / sqrt(3)
    
    # Calculate laminant realization
    volume_fraction = runif(1, volume_fraction_nominal-volume_fraction_variation, volume_fraction_nominal+volume_fraction_variation)
    volume_fraction = truncate_between(volume_fraction, 0, 1)
    
    # Form the constituents section
    constituents <- list("*CONSTITUENTS", "  NMATS=2", paste0("  ", parse_row(mat1)), paste0("  ", parse_row(mat2)))
    
    for (li in 1:length(layups)){
      clam <- eval(as.symbol(paste0("lam", li)))
      
      ## Form simple laminate
      n_ply <- 24
      cangles <- rnorm(length(clam), mean=clam, sd=angle_degree_variation)
      lam_def <- data.frame(LY=1:n_ply, MOD=2, THK=1/n_ply, ANG=cangles, ARCHID=1, VF=volume_fraction, "F"=1, M=2)
      laminate <- append(list("*LAMINATE", sprintf("  NLY=%d", n_ply)), as.list(paste0("   ", parse_df(lam_def))))
      
      ## Apply a fixed mechanical load
      mech <- data.frame(NPT=rep(2L, 6), TI=rep("0.,1.", 6L), MAG=c(sprintf("0.,%.2f", max_strain), rep("0.,0.", 5L)), 
                         MODE=c(1, rep(2, 5)))
      mload <- append(list("*MECH", "  LOP=99"), as.list(parse_df(mech, pre_space = "  ")))
      
      ## Static failure analysis 
      criteria1 <- "   CRIT=%d X11=%f COMPR=SAM ACTION=1"
      criteria2 <- "   CRIT=%d X11=%f X22=%f X33=%f X23=%f X13=%f X12=%f COMPR=SAM ACTION=1"
      sfa <- list("*FAILURE_SUBCELL", " NMAT=2", 
                  "  MAT=1 NCRIT=1", sprintf(criteria1, 1, UTS_fiber), 
                  "  MAT=2 NCRIT=1", sprintf(criteria2, 6, UTS_matrix, UTS_matrix, UTS_matrix, UTS_matrix_shear, UTS_matrix_shear, UTS_matrix_shear))
      
      ## Solver 
      solver_settings <- list(METHOD=1, NPT=2, TI="0., 1.", STP=0.005, ITMAX=20, ERR=1e-4)
      solver <- list("*SOLVER", paste0("  ", parse_row(solver_settings)))
      
      ## Output section
      out <- list("*PRINT", "  NPL=3")
      
      ## Get the stress strain output as well 
      ss_out <- list("*XYPLOT", "  FREQ=1", "  LAMINATE=1", 
                     sprintf("  NAME=graphite_epoxy_ss%09d X=1 Y=10", cc), 
                     "  MACRO=0", "  MICRO=0")
      
      ## Generate test mac files
      props = list("#INPUT PROPERTIES",
                   sprintf("# Layup: %s", layups[li]),
                   sprintf("# E_fiber: %f (%d \u00b1 %d%%)", E_fiber, E_fiber_nominal, E_fiber_percent_variation),
                   sprintf("# UTS_fiber: %f (%d \u00b1 %d%%)", UTS_fiber, UTS_fiber_nominal, UTS_fiber_percent_variation),
                   sprintf("# nu_fiber: %#.3f (%#0.3f \u00b1 %d%%)", nu_fiber, nu_fiber_nominal, nu_fiber_percent_variation),
                   sprintf("# G_fiber: %f", G_fiber),
                   sprintf("# E_matrix: %f (%d \u00b1 %d%%)", E_matrix, E_matrix_nominal, E_matrix_percent_variation),
                   sprintf("# UTS_matrix: %f (%d \u00b1 %d%%)", UTS_matrix, UTS_matrix_nominal, UTS_matrix_percent_variation),
                   sprintf("# nu_matrix: %#.3f (%#0.3f \u00b1 %d%%)", nu_matrix, nu_matrix_nominal, nu_matrix_percent_variation),
                   sprintf("# G_matrix: %f", G_matrix),
                   sprintf("# Volume_fraction: %f (%f \u00b1 %#.3f)", volume_fraction, volume_fraction_nominal, volume_fraction_variation),
                   "#"
      )
      mac_file <- c(props, constituents, laminate, mload, sfa, solver, out, ss_out, "*END")
      readr::write_lines(mac_file, path = sprintf(paste0(folder,"graphite_epoxy_%09d.mac"), cc))
      cc <- cc + 1L
    }
  }
}

# This function takes the mean and percent variation and 
# returns a specific random value from a normal or uniform distribution
property_realization <- function(mean, percent_variation, uniform=TRUE, n=1){
  if (uniform){
    # Find low/high in case mean is negative
    low = min(mean-mean*percent_variation/100, mean+mean*percent_variation/100)
    high = max(mean-mean*percent_variation/100, mean+mean*percent_variation/100)
    runif(n, low, high)
  } else {
    rnorm(n, mean, abs(mean*percent_variation/100))
  }
}

generateMACold <- function(constituents, analysis_type, loading, failure){
  

  ## Graphite fiber 
  mat1 <- list(M=1, CMOD=6, MATID="U", MATDB=1)
  
  E_A <- rnorm(n = 1, mean = 388200, sd = 388200*0.03)
  fiber_str <- rnorm(n = 1, mean = 3500, sd = 3500*0.05)
  
  E_Am <- rnorm(n = 1, mean = 3450, sd = 3450*0.03)
  matrix_str <- rnorm(n = 1, mean = 80, sd = 80*0.05)
  
  G_A <- rnorm(n = 1, mean = 14900, sd= 14900*0.03)
  V <- E_A/2*G_A -1 
  mat1 <- append(mat1, list(EL=c(E_A,E_A,V,V,G_A,-0.68,9.74))) # moduli in MPa
  
  E_Am <- rnorm(n = 1, mean = 3450, sd = 3450*0.03)
  G_Am <- rnorm(n = 1, mean = 1280, sd = 1280*0.03)
  Vm <- E_Am/2*G_Am -1 
  mat2 <- list(M=2, CMOD=6, MATID="U", MATDB=1)
  mat2 <- append(mat2, list(EL=c(E_Am,E_Am,Vm,Vm,G_Am,45,45))) # moduli in MPa
  
  #   ## Material 1 (Fiber)
#   mat1$E_A <- rnorm(n = 1, mean = 388200, sd = 388200*0.03)
#   mat1$G_A <- rnorm(n = 1, mean = 14900, sd= 14900*0.03)
#   sc1$X11 <- rnorm(n = 1, mean = 3500, sd = 3500*0.05) # in Mega pascals
#   sc1$XC11 <- sc1$X11
#   
#   ## Material 2 (Matrix)
#   mat2$E_A <- rnorm(n = 1, mean = 3450, sd = 3450*0.03)
#   mat2$G_A <- rnorm(n = 1, mean = 1280, sd = 1280*0.03)
#   sc2$X11 <- rnorm(n = 1, mean = 80, sd = 80*0.05)
  
  ## Epoxy Matrix 
  # mat2 <- list(M=2, CMOD=6, MATID="U", MATDB=1)
  # mat2 <- append(mat2, list(EL=c(3450,3450,0.35,0.35,1278,-45.E-6,45E-6))) # moduli in MPa
  #   
  ## Form the constituents section
  constituents <- list("*CONSTITUENTS", "  NMATS=2", paste0("  ", parse_row(mat1)), paste0("  ", parse_row(mat2)))
  
  for (li in 1:length(layups)){
    clam <- eval(as.symbol(paste0("lam", li)))
    
    ## Form simple laminate
    n_ply <- 24
    lam_vf <- rnorm(n_ply, mean = 0.60, sd = 0.60*0.03)
    cangles <- rnorm(length(clam), mean = clam, sd = abs(clam*0.01))
    lam_def <- data.frame(LY=1:n_ply, MOD=2, THK=1/n_ply, ANG=cangles, ARCHID=1, VF=lam_vf, "F"=1, M=2)
    laminate <- append(list("*LAMINATE", sprintf("  NLY=%d", n_ply)), as.list(paste0("   ", parse_df(lam_def))))
    
    ## Apply a fixed mechanical load
    mech <- data.frame(NPT=rep(2L, 6), TI=rep("0.,1.", 6L), MAG=c("0.,0.05", rep("0.,0.", 5L)), 
                       MODE=c(1, rep(2, 5)))
    mload <- append(list("*MECH", "  LOP=99"), as.list(parse_df(mech, pre_space = "  ")))
    
    ## Static failure analysis 
    criteria <- "   CRIT=%d X11=%f COMPR=SAM ACTION=1"
    sfa <- list("*FAILURE_SUBCELL", " NMAT=2", 
                "  MAT=1 NCRIT=1", sprintf(criteria, 1, fiber_str), 
                "  MAT=2 NCRIT=1", sprintf(criteria, 1, matrix_str))
    
    ## Solver 
    solver_settings <- list(METHOD=1, NPT=2, TI="0., 1.", STP=0.001, ITMAX=20, ERR=1e-4)
    solver <- list("*SOLVER", paste0("  ", parse_row(solver_settings)))
    
    ## Output section
    out <- list("*PRINT", "  NPL=3")
    
    ## Get the stress strain output as well 
    ss_out <- list("*XYPLOT", "  FREQ=1", "  LAMINATE=1", 
                   sprintf("  NAME=graphite_epoxy_ss%d X=1 Y=10", cc), 
                   "  MACRO=0", "  MICRO=0")
  
  mat_str <- constituent_lines
  
  }
  # readr::write_lines(x = mat_str, path = "test.MAC")
  # NMATS=2
  # # -- Carbon fiber
  # M=1 CMOD=6 MATID=U MATDB=1 &
  #   EL=230995.65703,238326.98225,0.20000,0.20000,97500.00000,0.0E-6,0.0E-6
  # # -- Epoxy matrix
  # M=2 CMOD=6 MATID=U MATDB=1 &
  #   EL=3454.93265,3454.93265,0.35403,0.35403,1275.79385,0.0E-6,0.0E-6
}

#' parse_constituents
#' @description Validates a data.frame of constituents for simulation with MAC/GMC. The data.frame is 
#' returned invisibly, and a preview of what the section would look like when produced with 
#' \code{\link{generateMAC}} is shown.  
#' @param constituents data.frame of constituent materials. See details. 
#' @details Must have at least 11 columns, plus optionally a name column.  
#' @references MAC/GMC Reference Manual v4.0: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20030015411.pdf
#' @export 
validate_constituents <- function(constituents){
  ## Validate data.frame
  stopifnot(is.data.frame(constituents))
  const_props <- c("M","CMOD","MATID","MATDB","E_A","E_T","V_A","V_T","G_A","ALPHA_A","ALPHA_T")
  stopifnot(all(toupper(colnames(constituents)) %in% const_props))
  
  ## Insert other checks for e.g. material property validity...
  
  ## Preview parsed MAC/GMC section
  parse_constituents(constituents, preview=TRUE)
  
  ## Return data.frame invisibly 
  invisible(constituents)
}

parse_constituents <- function(constituents, preview=FALSE){
  colnames(constituents) <- toupper(colnames(constituents))
  
  ## Parse non-elastic material properties 
  non_props <- parse_df(constituents[,c("M", "CMOD", "MATID", "MATDB")])
  
  ## Parse elastic properties
  el_props <- c("E_A","E_T","V_A","V_T","G_A","ALPHA_A","ALPHA_T")
  el_strs <- apply(constituents[,el_props], 1, function(el) paste0(el,collapse=", "))
  
  ## Combine 
  all_props <- mapply(function(np, p){ paste(np, p) }, non_props, as.list(paste0("EL=", el_strs)))
  
  constituent_lines <- list("*CONSTITUENTS", NMATS = nrow(constituents)) 
  
}
  
## Example output
# *CONSTITUENTS
#   NMATS=2
#   M=1 CMOD=6 MATID=U MATDB=1 EL=395552.88521, 395552.88521, 2849060785.65591, 2849060785.65591, 14405.46078, -0.68000, 9.74000
#   M=2 CMOD=6 MATID=U MATDB=1 EL=3476.48057, 3476.48057, 2274470.46133, 2274470.46133, 1308.49082, 45.00000, 45.00000
# *LAMINATE
#   NLY=24
#    LY=1 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.56449 F=1 M=2
#    LY=2 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.61081 F=1 M=2
#    LY=3 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.58791 F=1 M=2
#    LY=4 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.61394 F=1 M=2
#    LY=5 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.64251 F=1 M=2
#    LY=6 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.59719 F=1 M=2
#    LY=7 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.59674 F=1 M=2
#    LY=8 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.60497 F=1 M=2
#    LY=9 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.56399 F=1 M=2
#    LY=10 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.60731 F=1 M=2
#    LY=11 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.59035 F=1 M=2
#    LY=12 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.62355 F=1 M=2
#    LY=13 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.63239 F=1 M=2
#    LY=14 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.60966 F=1 M=2
#    LY=15 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.60743 F=1 M=2
#    LY=16 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.60002 F=1 M=2
#    LY=17 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.59181 F=1 M=2
#    LY=18 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.62730 F=1 M=2
#    LY=19 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.59576 F=1 M=2
#    LY=20 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.58345 F=1 M=2
#    LY=21 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.59152 F=1 M=2
#    LY=22 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.63803 F=1 M=2
#    LY=23 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.62079 F=1 M=2
#    LY=24 MOD=2 THK=0.04167 ANG=0 ARCHID=1 VF=0.60768 F=1 M=2
# *MECH
#   LOP=99
#   NPT=2 TI=0.,1. MAG=0.,0.05 MODE=1
#   NPT=2 TI=0.,1. MAG=0.,0. MODE=2
#   NPT=2 TI=0.,1. MAG=0.,0. MODE=2
#   NPT=2 TI=0.,1. MAG=0.,0. MODE=2
#   NPT=2 TI=0.,1. MAG=0.,0. MODE=2
#   NPT=2 TI=0.,1. MAG=0.,0. MODE=2
# *FAILURE_SUBCELL
#  NMAT=2
#   MAT=1 NCRIT=1
#    CRIT=1 X11=3545.645689 COMPR=SAM ACTION=1
#   MAT=2 NCRIT=1
#    CRIT=1 X11=78.935982 COMPR=SAM ACTION=1
# *SOLVER
#   METHOD=1 NPT=2 TI=0., 1. STP=0.00100 ITMAX=20 ERR=0.00010
# *PRINT
#   NPL=3
# *XYPLOT
#   FREQ=1
#   LAMINATE=1
#   NAME=graphite_epoxy_ss1 X=1 Y=10
#   MACRO=0
#   MICRO=0
# *END

  
  
  