zinck.filter.interaction <- function(X, X_tilde, Y, Z, model = "Random Forest", 
                                     fdr = 0.1, offset = 1, seed = NULL, nlambda = 500,
                                     ntrees = 1000, tune_mtry = FALSE,  mtry = NULL, metric = NULL,
                                     rftuning = FALSE, conservative = FALSE, 
                                     normalization = NULL, heredity = "strong") {
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
  
  # ========== MODEL FITTING ==========
  if (model == "glmnet") {
    # Use glmnet for lasso regression
    W_stats <- fit_glmnet_interaction(X, X_tilde, X_Z, X_tilde_Z, Y, Z,
                                      nlambda = nlambda, seed = seed)
    W_main <- W_stats$W_main
    W_int <- W_stats$W_int
  } else if (model == "Random Forest"){
    if (rftuning == TRUE) {
      W_stats <- fit_rf_interaction(X, X_tilde, X_Z, X_tilde_Z, Y, Z,
                                    ntrees = ntrees, tune_mtry = tune_mtry,
                                    mtry = mtry, metric = metric, seed = seed)
      W_main <- W_stats$W_main
      W_int <- W_stats$W_int
    } else { # Use built-in knockoff statistics functions
      if (!is.null(seed)){set.seed(seed)}
      W_main <- stat.random_forest(X, X_tilde, Y)
      W_int <- stat.random_forest(X_Z, X_tilde_Z, Y)
      names(W_main) <- colnames(X)
      names(W_int) <- colnames(X)
    }
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
