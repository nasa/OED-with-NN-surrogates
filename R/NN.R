#' MLP 
#' @description Creates a multilayer perceptron using \pkg{keras}'s functional API. 
#' @param layers named-list of layer parameters, whose names correspond to layer names and value to layer parameters. 
#' @details The \pkg{keras} package simplifies building neural networks in R using the 
#' \pkg{magrittr} pipe-compatible API to functionally compose layers. But individual layers 
#' are linked at the time of evaluation, making it cumbersome to e.g. programmatically change 
#' parameters at the individual layer level. This function uses the magrittr pipe ('%>%') 
#' to turn a simple named-list of layers + parameters into a functional sequence of activation layers. 
#' @return a unary function right-folded using the magrittr pipe whose evaluation yields a keras model. 
#' @import magrittr reticulate tensorflow keras
#' @examples
#' 
#' ## Make simple regression model
#' X <- matrix(runif(100))
#' Y <- (0.5 * X + 2) + rnorm(100, sd=0.2)
#' layers <- list(
#'  "dense" = list(units = 1, activation = 'linear')
#' )
#' 
#' ## Programmatically change the parameters 
#' for (i in seq(10)){
#'   layers[[1]]$units <- i
#'   model <- MLP(layers)
#'   ## ... 
#' }
#' 
MLP <- function(layers){
  mlp <- lapply(1:length(layers), function(i){
    f <- eval(as.symbol(paste0("layer_", names(layers)[[i]])))
    (do.call(f, layers[[i]]))
  })
  # function(x) {magrittr::freduce(x, mlp) }
  . %>% { magrittr::freduce(., mlp) }

  # model <- keras::keras_model_sequential()
  # for(i in 1:length(layers)) {
  #   layer_f <- eval(as.symbol(paste0("layer_", names(layers)[[i]])))
  #   if (i == 1L){ layers[[i]]$input_shape = n_properties }
  #   if (i == length(layers)){ layers[[i]]$units = n_classes }
  #   model$add(do.call(layer_f, layers[[i]]))
  # }
}

#' Partition 
#' @description Partitions a data set into training, testing, and validation splits
#' @param X data.table of input variables.
#' @param Y data.table of output variables. 
#' @return List with partitioned indices for the training, testing, and validation sets. 
partition <- function(X, Y, input_shape, splits = c(train = 0.80, validate = 0.10, test = 0.10)){
  if (!is(X, "data.table")){ stop("'MLP' expects 'X' to be a named data.table.") }
  # if (!is(Y, "integer")){ stop("'MLP' expects 'Y' to be an integer vector of classes.") }
  
  set.seed(1234) ## for reproducibility
  if (is(Y, "integer")){
    ## One-hot encoding assumes 0-based class labels
    Y <- match(Y, unique(Y)) - 1L
    
    ## Choose the subset to use for training
    class_labels <- unique(Y)
    n_classes <- length(class_labels)
    
    ## Split into training, testing, and validation, equally among classes
    class_splits <- lapply(class_labels, function(cid){
      class_idx <- which(Y == cid)
      class_len <- length(class_idx)
      class_partition <- cut(seq(class_len), breaks = class_len * cumsum(c(0, splits)), labels = names(splits))
      split(class_idx, sample(class_partition))
    })
  } else if (is(Y, "matrix")){
    class_partition <- cut(seq(nrow(X))/nrow(X), breaks = cumsum(c(0, splits)), labels = names(splits))
    class_splits <- split(seq(nrow(X)), sample(class_partition))
  }
  
  
  ## Make the training, validation, and testing sets
  if (is(Y, "integer")){
    train_idx <- as.vector(unlist(sapply(class_splits, function(spl){ spl[["train"]] })))
    test_idx <- as.vector(unlist(sapply(class_splits, function(spl){ spl[["test"]] })))
    validate_idx <- as.vector(unlist(sapply(class_splits, function(spl){ spl[["validate"]] })))
  } else if (is(Y, "matrix")){
    train_idx <- as.vector(class_splits[["train"]])
    test_idx <- as.vector(class_splits[["test"]])
    validate_idx <- as.vector(class_splits[["validate"]])
  }
  
  ## Useful variables 
  n_train <- length(train_idx)
  n_test <- length(test_idx)
  n_validate <- length(validate_idx)
  
  ## Separate training and testing data 
  x_train <- array_reshape(x = as.matrix(X[train_idx,]), dim = c(n_train, input_shape))
  x_test <- array_reshape(x = as.matrix(X[test_idx,]), dim = c(n_test, input_shape))
  x_val <- array_reshape(x = as.matrix(X[validate_idx,]), dim = c(n_validate, input_shape))
  
  if (is(Y, "integer")){
    ## One-hot encode Y 
    y_train <- to_categorical(Y[train_idx], num_classes = n_classes)
    y_test <- to_categorical(Y[test_idx], num_classes = n_classes)
    y_val <- to_categorical(Y[validate_idx], num_classes = n_classes)
  } else if (is(Y, "matrix")){
    ## Just subset
    y_train <- Y[train_idx,]
    y_test <- Y[test_idx,]
    y_val <- Y[validate_idx,]
  }
  
  ## Return the model
  return(list(train = list(X = x_train, Y = y_train, idx = train_idx), 
              val = list(X = x_val, Y = y_val, idx = validate_idx), 
              test = list(X = x_test, Y = y_test, idx = test_idx)))
}


