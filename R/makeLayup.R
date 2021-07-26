#' makeLayup 
#' @description generates a ply layup using standard notation. 
#' @param x layup orientation. See details. 
#' @export
#' @examples 
#' makeLayup("[0]_24")
#' makeLayup("[0/90]_12")
#' makeLayup("[0, 15, 75, 90]_3S")
makeLayup <- function(x){
  inner <- gsub(x, pattern = "\\[(.*)\\].*", replacement = "\\1")
  angles <- as.numeric(strsplit(inner, split = "[,/]")[[1]])
  outer <- gsub(x, pattern = "\\[.*\\]_(.*)", replacement = "\\1")
  n <- as.integer(gsub(outer, pattern = "(\\d+)S?", replacement = "\\1"))
  if (length(grep(pattern = "S", x = outer, ignore.case = TRUE)) > 0){
    rep(c(angles, rev(angles)), n)
  } else {
    rep(angles, n)
  }
}