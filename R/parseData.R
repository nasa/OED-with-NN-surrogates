## Simple function to retrieve subsets of the data, with the option of scaling to (zero mean, unit sd), 
## converting dimensions to numeric (such as categorical or ordinal features), and whether a data matrix 
## of data table is requested. 
#' @export
parseData <- function(x, dim = NULL, scale_col=TRUE, remove.NA=TRUE, all_numeric = FALSE, as_matrix=FALSE){
  if (!is(x, "data.table")) { "'parseData' expects 'x' to be a data.table object. "}
  if (!is(dim, "character")) { "'parseData' expects 'dim' to be a character vector. "}
  if (missing(dim)){ dim <- colnames(x) }
  tmp <- x[, ..dim]
  if (remove.NA){ tmp <- na.omit(tmp) }
  if (all_numeric){ tmp <- tmp[, lapply(.SD, dummy_code)] }
  if (scale_col){ tmp <- tmp[, lapply(.SD, function(x_i) if(is.numeric(x_i)) scale(x_i) else x_i)] }
  return(if(as_matrix){ as.matrix(tmp) } else { tmp })
}


#' @export
dummy_code <- function(x){ if (is.factor(x)){ as.integer(x) } else { as.numeric(x) } }