# stat_rf_interaction.R
# r package ranger is required!
## revise the built-in rf function in knockoff package (we need to conduct rf on the whole augmented design).
stat_rf_interaction <- function(X, X_tilde, X_Z, X_tilde_Z, Y, Z, seed = NULL) {
  p <- ncol(X)
  n <- nrow(X)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  swap <- rbinom(p, 1, 0.5)
  swap.M <- matrix(swap, nrow = n, ncol = p, byrow = TRUE)
  
  X_swap <- X * (1 - swap.M) + X_tilde * swap.M
  X_tilde_swap <- X * swap.M + X_tilde * (1 - swap.M)
  X_Z_swap <- X_Z * (1 - swap.M) + X_tilde_Z * swap.M
  X_tilde_Z_swap <- X_Z * swap.M + X_tilde_Z * (1 - swap.M)
  
  X_aug <- cbind(Z, X_swap, X_Z_swap, X_tilde_swap, X_tilde_Z_swap)
  
  df <- data.frame(y = Y, X_aug)
  rfFit <- ranger::ranger(y ~ ., data = df, importance = "impurity", 
                          write.forest = FALSE, num.trees = 1000)
  imp <- as.vector(rfFit$variable.importance)

  imp_X <- imp[2:(1 + p)]
  imp_X_Z <- imp[(2 + p):(1 + 2*p)]
  imp_X_tilde <- imp[(2 + 2*p):(1 + 3*p)]
  imp_X_tilde_Z <- imp[(2 + 3*p):(1 + 4*p)]
  
  W_main <- (abs(imp_X) - abs(imp_X_tilde)) * (1 - 2*swap)
  W_int <- (abs(imp_X_Z) - abs(imp_X_tilde_Z)) * (1 - 2*swap)
  
  names(W_main) <- colnames(X)
  names(W_int) <- colnames(X)
  
  return(list(W_main = W_main, W_int = W_int))
}
