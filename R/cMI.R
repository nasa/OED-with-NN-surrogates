#' Fast conditional Mutual Information
#' @description Computes the conditional mutual information between three random variables, i.e. I(X ; Y | Z)
#' Uses a fast range search to compute the neighbors in the max-norm. 
#' Relies on Kraskov's estimator
#' @import kdtools
#' @export 
cMI <- function(X, Y, Z, k = "suggest", alpha = 0.7, type = c("ps", "ksg")){
  N <- nrow(X)
  { d_x <- ncol(X); d_y <- ncol(Y); d_z <- ncol(Z) }
  { c_dx <- 2^d_x; c_dy <- 2^d_y; c_dx <- 2^d_z }
  if(missing(k) || k == "suggest"){ k <- ifelse((d_x+d_y+d_z)*log(N) >= N, 15L, as.integer((d_x+d_y+d_z)*log(N))) }
  
  if (missing(type) || type == "ps"){
    v_d <- function(d){ pi^(d/2)/gamma(d/2 + 1) }
    { d_xyz <- d_x+d_y+d_z; d_xz <- d_x+d_z; d_yz <- d_y+d_z }
    { c_xyz <- v_d(d_x+d_y+d_z); c_xz <- v_d(d_x+d_z); c_yz <- v_d(d_y+d_z); c_z <- v_d(d_z) }
    B <- (gamma(k)^2)/(gamma(k-alpha+1)*gamma(k+alpha-1))
    ai <- (1-alpha)
    n1 <- RANN::nn2(data = cbind(X, Y, Z), k = k+1L)$nn.dist[, k+1L]^(d_xyz*ai)
    n2 <- RANN::nn2(data = cbind(X, Z), k = k+1L)$nn.dist[, k+1L]^(d_xz*ai)
    n3 <- RANN::nn2(data = cbind(Y, Z), k = k+1L)$nn.dist[, k+1L]^(d_yz*ai)
    n4 <- RANN::nn2(data = Z, k = k+1L)$nn.dist[, k+1L]^(d_z*ai)
    return(sum((n1/n2)*(n3/n4)*B^2)/N)
  }

    
  if (type == "ksg"){
    nn <- rmi::nearest_neighbors(cbind(X, Y, Z), k = k)
    eps <- nn$nn_dist[, k+1] + sqrt(.Machine$double.eps)
    
    ## Compute the number of neighbors in each max-norm rectangle in the marginal spaces
    subspace_neighbors <- function(M){
      M_ <- matrix_to_tuples(M)
      M_idx <- kd_order(M_)
      kd_sort(M_, inplace=TRUE)
      sapply(1:N, function(i){ 
        mpt <- M_[M_idx[i],]
        nrow(kd_range_query(M_, mpt-(eps[i]/2), mpt+(eps[i]/2)))
      })
    }
    
    # { XZ <- matrix_to_tuples(cbind(X, Z));  YZ <- matrix_to_tuples(cbind(Y, Z)); Z_ <- matrix_to_tuples(Z) }
    # { xz_idx <- kd_order(XZ); yz_idx <- kd_order(YZ); z_idx <- kd_order(Z) }
    # kd_sort(XZ, inplace = TRUE)
    # kd_sort(YZ, inplace = TRUE)
    # kd_sort(Z_, inplace = TRUE)
    # n_xz <- sapply(1:N, function(i){ nrow(kd_range_query(XZ, XZ[xz_idx[i],]-(eps[i]/2), XZ[xz_idx[i],]+(eps[i]/2)))  })
    # n_yz <- sapply(1:N, function(i){ nrow(kd_range_query(YZ, YZ[yz_idx[i],]-(eps[i]/2), YZ[yz_idx[i],]+(eps[i]/2)))  })
    # n_z <- sapply(1:N, function(i){ nrow(kd_range_query(Z_, Z_[z_idx[i],]-(eps[i]/2), Z_[z_idx[i],]+(eps[i]/2)))  })
    
    ## Expected value of rectilinear box + correction term
    d <- d_x+d_y+d_z
    e_corr <- (d-1)*digamma(N)+digamma(k)-(d-1)/k
    
    ## Neighbors in the marginal spaces
    nx <- subspace_neighbors(X)
    ny <- subspace_neighbors(Y)
    nz <- subspace_neighbors(Z)
    
    ## I(X;Y;Z)
    I_xyz <- e_corr - (1/N)*sum(c(digamma(nx), digamma(ny), digamma(nz)))
    
    ## I(X;Z)
    I_xz <- e_corr - (1/N*(sum(c(digamma(nx), digamma(nz)))))
    
    ## I(X;Y)
    I_xy <- e_corr - (1/N*sum(c(digamma(nx), digamma(ny))))
    
    ## H(X), H(Z)
    h_x <- digamma(N)-digamma(k)+(d_x/N)*sum(log(eps))
    h_z <- digamma(N)-digamma(k)+(d_z/N)*sum(log(eps))
    
    ## I(X;Y|Z)
    return(I_xyz-(I_xz+I_xy)+h_x-h_z)
  }

  
  # 
  # I <- function()
  # 
  # 
  # ## Harmonic-type estimate based on: "Partial Mutual Information for Coupling Analysis of Multivariate Time Series"
  # h_xz <- sapply(n_xz+1, function(nxz) -sum(1/1:nxz))
  # h_yz <- sapply(n_yz+1, function(nyz) -sum(1/1:nyz))
  # h_z <- sapply(n_z+1, function(nz) -sum(1/1:nz))
  # cmi_est <- mean(h_xz + h_yz - h_z) + sum(1/1:k)
  
  ## KSG estimator of CMI using RMI package 
  # I_xyz <- rmi::knn_mi(cbind(X, Y, Z), splits = c(1,1,1), options = list(method = "KSG1", k = k))
  # I_xz <- rmi::knn_mi(cbind(X, Z), splits = c(ncol(X), ncol(Z)), options = list(method = "KSG1", k = k))
  # I_xy <- rmi::knn_mi(cbind(X, Y), splits = c(ncol(X), ncol(Y)), options = list(method = "KSG1", k = k))
  # h_x <- digamma(N) - digamma(k) + d_x/N * sum(2*eps)
  # h_z <- digamma(N) - digamma(k) + d_z/N * sum(2*eps)
  # return(I_xyz -(I_xz + I_xy) + h_x - h_z)
  
  # 
  # mean(-sum(1/(n_xz + 1)) + -sum(1/(n_yz + 1)) - -sum(1/(n_z + 1)))
  # (1/N) * sum(digamma(N) + log(n_xz+1) + log(n_yz+1) + log(n_z+1))
}