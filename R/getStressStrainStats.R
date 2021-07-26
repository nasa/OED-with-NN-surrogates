#' getStressStrainStats
#' 
#' Attempts to detect various statistics along the stress-strain curve numerically.  
#' @export
getStressStrainStats <- function(strain, stress, uts_linearity_thresh = 0.98){
  n <- length(stress)
  
  ## Test linearity along every point in the ordered curve
  pearson_cc <- sapply(2:length(stress), function(i) cor(stress[1:i], strain[1:i]))
  pl_idx <- max(which(abs(pearson_cc - 1.0) < sqrt(.Machine$double.eps))) + 1 
  linear_idx <- which(pearson_cc > uts_linearity_thresh)
  #uts_idx <- which.max(stress[1:max(linear_idx)+1L])
  uts_idx <- which.max(stress)
  ## Return a named vector 
  res <- list(ult_tensile_stress=stress[uts_idx], 
              ult_tensile_strain=strain[uts_idx],
              proportionality_limit_stress=stress[pl_idx],
              proportionality_limit_strain=strain[pl_idx],
              elastic_modulus=stress[pl_idx]/strain[pl_idx])
  return(unlist(res))
}
