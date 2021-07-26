#include <Rcpp.h>
using namespace Rcpp;


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
layups <- c("[0]_24", 
            "[90]_24", 
            "[0/90]_12", 
            "[0/90]_6S", 
            "[15, -15]_12", 
            "[30, -30]_12",
            "[45, -45]_12",
            "[60, -60]_12", 
            "[75, -75]_12",
            "[15, -15]_6S",
            "[30, -30]_6S",
            "[45, -45]_6S",
            "[60, -60]_6S",
            "[75, -75]_6S",
            "[0, -45, 45, 90]_3S", 
            "[0, -60, 60]_4S", 
            "[0, 90, -30, 30, -60, 60]_2S", 
            "[0, 45, 90, -45]_6", 
            "[45, -45, 90, 0]_3S") 
for (li in 1:length(layups)){ assign(paste0("lam", li), value = makeLayup(layups[[li]])) }
lam_layups <- sapply(1:length(layups), function(i) eval(as.symbol(sprintf("lam%d", i))))
## Normalize closure
normalizer <- function(x, a = -1, b = 1, rng=NULL){ 
  # if (!is.null(opts) list("ab_range"=list(a=0, b=1), =list(center=TRUE, scale=TRUE))
  if (is.vector(x)){
    { min_x <- min(x); max_x <- max(x) }
    normalize <- function(x, invert=FALSE){ 
      if (missing(invert) || !invert){ return(((b - a)*(x - min_x)) / (max_x - min_x) + a ) }
      else { return(((x - a)/(b - a))*(max_x - min_x) + min_x) }
    }
    return(normalize)
  }
  else if (is.matrix(x)){
    d_x <- ncol(x)
    if (missing(rng) || is.null(rng)){ 
      col_idx <- as.list(seq(ncol(x)))
      x_rng <- apply(x, 2, range)
    } else { 
      if (sum(rng) != d_x) { stop("'rng' should sum to the number fo columns in x.")}
      start_idx <- cumsum(rng) - rng + 1
      col_idx <- lapply(1:length(rng), function(i) start_idx[i]:(start_idx[i]+rng[i]-1L))
      x_rng <- sapply(1:length(col_idx), function(i) { range(as.vector(unlist(x[, col_idx[[i]]]))) })  
    }
    normalize <- function(x, invert=FALSE){ 
      if (missing(invert) || !invert){ 
        x_norm <- matrix(0, ncol = ncol(x), nrow = nrow(x))
        for (i in 1:length(col_idx)){
          idx <- col_idx[[i]]
          x_tmp <- as.vector(unlist(x[, idx]))
          x_norm[, idx] <- ((b - a)*(x_tmp - x_rng[1, i])) / (x_rng[2, i] - x_rng[1, i]) + a 
        }
        return(x_norm) 
      }
      else {
        x_unnorm <- matrix(0, ncol = ncol(x), nrow = nrow(x))
        for (i in 1:length(col_idx)){
          idx <- col_idx[[i]]
          x_tmp <- as.vector(unlist(x[, idx]))
          x_unnorm[, idx] <- (((x_tmp - a)/(b - a))*(x_rng[2, i] - x_rng[1, i]) + x_rng[1, i]) 
        }
        return(x_unnorm)
      }
    }
  }
  return(normalize)
}
layup_normalizer <- normalizer(1:length(layups))


load("~/res2.rdata")
hist(unlist(lapply(1:23, function(i) { res[[i]]$phase1.d[,7] })))
library("animate")

load("~/res4.rdata")
cc <- 1L
i <- 1L
animation::saveGIF({
  freq <- NULL
  
  while( i < 23) {
    freq <- c(freq, round(layup_normalizer(res[[i]]$phase1.d[, 8], invert = TRUE)))
    freq_ord <- ordered(layups[freq], layups[c(1:17, 19, 18)])
    res_tbl <- table(freq_ord)
    barplot(res_tbl, main = "Optimal Design Laminate frequencies", las=2, cex.names = 0.6)
    i <- i + 1L
    if (i == 23 && cc < 3){
      i <- 1
      cc <- cc + 1
    }
  }
}, interval = 0.25)
lapply(1:23, function(i) { res[[1]]$phase1.d[,1] })




res_table <- table(layups[unlist(ods)])
x <- barplot(res_table, xaxt = "n", main = "Optimal Design Laminate frequencies")
# barplot(res_table, las=2, cex.names = 0.6)
text(cex=1, x=x-.25, y=-1.25, labs=labels(res_table)[[1]], xpd=TRUE, srt=45)



p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity")
p




### BAR PLOT 

load("~/res3.rdata")
ods <- lapply(1:25, function(i) {
  round(layup_normalizer(res[[i]]$phase1.d[, 6], invert = TRUE))
})

res_table <- table(layups[unlist(ods)])
library("ggplot2")
df <- as.data.frame(res_table)
colnames(df) <- c("Layup", "Frequency")
df$Layup <- ordered(df$Layup, layups[c(1:17, 19, 18)])

layup_colors <- viridis::viridis(length(df$Layup)) # rainbow(length(layups))

p <- ggplot(data = df, aes(x=Layup, y=Frequency, fill = Layup)) +
  # geom_text(aes(label=len), vjust=-0.3, size=3.5) + 
  geom_bar(stat="identity", color = layup_colors) + 
  ggtitle("Optimal Laminate Frequency (-fiber/matrix strengths)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) + 
  scale_fill_manual(values = layup_colors)
  # theme_minimal()


#########################


layups
p + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
df


res1_tbl <- table(layups[as.integer(res$phase2.d)])
barplot(res1_tbl, las=2, cex.names = 0.6)
# text(cex=1, x=1:length(res2_tbl), y=-1.25, labels=names(res2_tbl), xpd=TRUE, srt=90)
layups[which(!layups %in% names(res1_tbl))]

res2_tbl <- table(layups[as.integer(res2$phase2.d)])
barplot(res2_tbl, las=2, cex.names = 0.6)
# text(cex=1, x=1:length(res2_tbl), y=-1.25, labels=names(res2_tbl), xpd=TRUE, srt=90)
layups[which(!layups %in% names(res2_tbl))]


lapply(1:21, function(i) {
  table(round(layup_normalizer(res[[i]]$phase1.d[, 6], invert = TRUE)))
})

*/
