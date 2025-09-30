zinck.filter.interaction <- function(X, X_tilde, Y, Z, model = "Random Forest", 
                                     fdr = 0.1, offset = 1, seed = NULL, 
                                     ntrees = 1000, tune_mtry = FALSE,  mtry = NULL, metric = NULL,
                                     rftuning = FALSE, conservative = FALSE, 
                                     normalization = NULL, heredity = c("separate", "strong")) {
  if (!is.null(normalization)) {
    if (normalization == "prop"){
      X <- convert_to_proportions(X)
      X_tilde <- convert_to_proportions(X_tilde)
    } else if (normalization == "log"){
      X <- log_normalize(X)
      X_tilde <- log_normalize(X_tilde)
    } else if (normalization == "CLR"){
      X <- CLR_normalize(X)
      X_tilde <- CLR_normalize(X_tilde)
    }
  }
  
  X_Z <- X * Z  # Element-wise multiplication with Z (exposure)
  X_tilde_Z <- X_tilde * Z  # Knockoff interaction terms
  
  # Construct augmented design matrix: {Z, X, X^(Z), X_tilde, X_tilde^(Z)}
  X_aug <- cbind(Z, X, X_Z, X_tilde, X_tilde_Z)
  
  colnames(X_aug) <- c("Z", 
                       colnames(X), 
                       paste0(colnames(X), "_Z"), 
                       paste0(colnames(X_tilde), "_tilde"), 
                       paste0(colnames(X_tilde), "_Z_tilde"))
  
  if (model == "Random Forest") {
    if (rftuning == TRUE){
      if (is.null(mtry)) { # Tune mtry if not provided
        if (tune_mtry) {
          if (is.factor(Y) || length(unique(Y)) == 2) { # Binary response
            set.seed(seed)
            bestmtry <- tuneRF(X_aug, as.factor(Y), stepFactor = 1.5, improve = 1e-5, ntree = ntrees, trace = FALSE)
            mtry <- bestmtry[as.numeric(which.min(bestmtry[, "OOBError"])), 1]
        } else{
          set.seed(seed)
          bestmtry <- tuneRF(X_aug,Y, stepFactor = 1.5, improve = 1e-5, ntree = ntrees, trace = FALSE)
          mtry <- bestmtry[as.numeric(which.min(bestmtry[, "OOBError"])), 1]
      }}
      else {
        mtry <- floor(sqrt(ncol(X_aug))) # Default mtry
      }
    }
    if (is.factor(Y) || length(unique(Y)) == 2) { # Binary response
      set.seed(seed)
      model_rf <- randomForest(X_aug, as.factor(Y), ntree = ntrees, mtry = mtry, importance = TRUE)
      if (metric == "Accuracy"){
        cf <- randomForest::importance(model_rf)[, 1]
    } else if (metric == "Gini"){
      cf <- randomForest::importance(model_rf)[,3]
    }
    } else { # continuous response
      set.seed(seed)
      model_rf <- randomForest(X_aug, Y, ntree = ntrees, mtry = mtry, importance = TRUE)
      if (metric == "Accuracy"){ 
        cf <- randomForest::importance(model_rf)[, 1]
      } else if (metric == "Gini"){
        cf <- randomForest::importance(model_rf)[,3]
      }
    }
      
    imp_X <- cf[2:(1 + ncol(X))]                    # Main effects
    imp_X_Z <- cf[(2 + ncol(X)):(1 + 2*ncol(X))]   # Interaction effects
    imp_X_tilde <- cf[(2 + 2*ncol(X)):(1 + 3*ncol(X))]  # Knockoff main
    imp_X_tilde_Z <- cf[(2 + 3*ncol(X)):(1 + 4*ncol(X))] # Knockoff interaction
    
    W_main <- abs(imp_X) - abs(imp_X_tilde)
    W_int <- abs(imp_X_Z) - abs(imp_X_tilde_Z)
    
    names(W_main) <- colnames(X)
    names(W_int) <- colnames(X)
    
    } else if (rftuning == FALSE) {
      # Use built-in knockoff statistics functions
      set.seed(seed)
      # Main effects using built-in function
      if (is.factor(Y) || length(unique(Y)) == 2) { # binary
        W_main <- stat.random_forest(X, X_tilde, as.factor(Y))
      } else { # continous
        W_main <- stat.random_forest(X, X_tilde, Y)
      }
      names(W_main) <- colnames(X)
      
      # Interaction effects using built-in function
      if (is.factor(Y) || length(unique(Y)) == 2) {
        W_int <- stat.random_forest(X_Z, X_tilde_Z, as.factor(Y))
      } else {
        W_int <- stat.random_forest(X_Z, X_tilde_Z, Y)
      }
      names(W_int) <- colnames(X)
    }
    
    # Apply selection procedure based on heredity constraint
    if (heredity == "separate") {
      # Option A: Separate filters
      if (conservative == FALSE) {
        T_main <- knockoff.threshold(W_main, fdr = fdr, offset = offset)
        T_int <- knockoff.threshold(W_int, fdr = fdr, offset = offset)
      } else {
        T_main <- (1 - fdr) * ko.sel(W_main, print = FALSE, method = "gaps")$threshold
        T_int <- (1 - fdr) * ko.sel(W_int, print = FALSE, method = "gaps")$threshold
      }
      
      selected_main <- sort(which(W_main >= T_main))
      selected_int <- sort(which(W_int >= T_int))
      
      out <- list(
        selected_main = selected_main,
        selected_interaction = selected_int,
        W_main = W_main,
        W_int = W_int,
        T_main = T_main,
        T_int = T_int,
        heredity = "separate"
      )
      
    } else if (heredity == "strong") {
      # Option B: Two-layer hierarchical knockoffs for strong heredity
      
      # Layer 1: Main effects
      if (conservative == FALSE) {
        T_main <- knockoff.threshold(W_main, fdr = fdr, offset = offset)
      } else {
        T_main <- (1 - fdr) * ko.sel(W_main, print = FALSE, method = "gaps")$threshold
      }
      
      selected_main <- sort(which(W_main >= T_main))
      
      # Layer 2: Interaction effects restricted to selected main effects
      if (length(selected_main) > 0) {
        W_int_restricted <- W_int[selected_main]
        
        if (conservative == FALSE) {
          T_int <- knockoff.threshold(W_int_restricted, fdr = fdr, offset = offset)
        } else {
          T_int <- (1 - fdr) * ko.sel(W_int_restricted, print = FALSE, method = "gaps")$threshold
        }
        
        selected_int <- selected_main[which(W_int_restricted >= T_int)]
      } else {
        selected_int <- integer(0)
        T_int <- Inf
      }
      
      out <- list(
        selected_main = selected_main,
        selected_interaction = selected_int,
        W_main = W_main,
        W_int = W_int,
        T_main = T_main,
        T_int = T_int,
        heredity = "strong"
      )
    }
  return(out)
  }
}
