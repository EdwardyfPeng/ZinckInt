## Using random forest to compute W
# r package randomForest is required!

fit_rf_interaction <- function(X, X_tilde, X_Z, X_tilde_Z, Y, Z,
                                ntrees = 1000, tune_mtry = FALSE, mtry = NULL,
                                metric = NULL, seed = NULL) {
  # Construct augmented design matrix
  X_aug <- cbind(Z, X, X_Z, X_tilde, X_tilde_Z)
  colnames(X_aug) <- c("Z",
                       colnames(X),
                       paste0(colnames(X), "_Z"),
                       paste0(colnames(X_tilde), "_tilde"),
                       paste0(colnames(X_tilde), "_Z_tilde"))
  # Tune or set mtry
  if (is.null(mtry)) {
    if (tune_mtry) {
      if (is.factor(Y) || length(unique(Y)) == 2) {
        set.seed(seed)
        bestmtry <- randomForest::tuneRF(X_aug, as.factor(Y), stepFactor = 1.5, 
                                         improve = 1e-5, ntree = ntrees, trace = FALSE)
        mtry <- bestmtry[as.numeric(which.min(bestmtry[, "OOBError"])), 1]
      } else {
        set.seed(seed)
        bestmtry <- randomForest::tuneRF(X_aug, Y, stepFactor = 1.5, 
                            improve = 1e-5, ntree = ntrees, trace = FALSE)
        mtry <- bestmtry[as.numeric(which.min(bestmtry[, "OOBError"])), 1]
      }
    } else { # otherwise use default mtry
      mtry <- floor(sqrt(ncol(X_aug)))
    }
  }
  
  # Fit random forest
  if (is.factor(Y) || length(unique(Y)) == 2) { # Binary response
    set.seed(seed)
    model_rf <- randomForest::randomForest(X_aug, as.factor(Y), ntree = ntrees, 
                                mtry = mtry, importance = TRUE)
    if (metric == "Accuracy"){
      cf <- randomForest::importance(model_rf)[, 1]
    } else if (metric == "Gini"){
      cf <- randomForest::importance(model_rf)[, 3]
    } else {
      cf <- randomForest::importance(model_rf)[, 3]  # default to Gini
    }
  } else { # Continuous response
    set.seed(seed)
    model_rf <- randomForest::randomForest(X_aug, Y, ntree = ntrees, 
                                mtry = mtry, importance = TRUE)
    if (metric == "Accuracy"){
      cf <- randomForest::importance(model_rf)[, 1]
    } else {
      cf <- randomForest::importance(model_rf)[, 1]  # default to Permutation Importance
    }
  }
  
  # Extract importance
  imp_X <- cf[2:(1 + ncol(X))]
  imp_X_Z <- cf[(2 + ncol(X)):(1 + 2*ncol(X))]          
  imp_X_tilde <- cf[(2 + 2*ncol(X)):(1 + 3*ncol(X))]
  imp_X_tilde_Z <- cf[(2 + 3*ncol(X)):(1 + 4*ncol(X))]
  
  # Compute W statistics
  W_main <- abs(imp_X) - abs(imp_X_tilde)
  W_int <- abs(imp_X_Z) - abs(imp_X_tilde_Z)
  
  names(W_main) <- colnames(X)
  names(W_int) <- colnames(X)
  
  return(list(W_main = W_main, W_int = W_int))
}
