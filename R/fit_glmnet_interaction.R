## Using Lasso regression to compute W
# r package "glmnet" is required!

fit_glmnet_interaction <- function(X, X_tilde, X_Z, X_tilde_Z, Y, Z, 
                                   nlambda = 500, seed = NULL) {
  p <- ncol(X)
  n <- nrow(X)
  
  if (!is.null(seed)){
    set.seed(seed)
  }

  # main effect terms and interaction terms should share the same swapping vector. 
  swap_main <- rbinom(p, 1, 0.5)
  swap_main.M <- matrix(swap_main, nrow=n, ncol=p, byrow=TRUE)
  X_swap <- X * (1 - swap_main.M) + X_tilde * swap_main.M
  X_tilde_swap <- X * swap_main.M + X_tilde * (1 - swap_main.M)
  
  X_Z_swap <- X_Z * (1 - swap_main.M) + X_tilde_Z * swap_main.M
  X_tilde_Z_swap <- X_Z * swap_main.M + X_tilde_Z * (1 - swap_main.M)
  
  X_aug <- cbind(Z, X_swap, X_Z_swap, X_tilde_swap, X_tilde_Z_swap)
  X_aug_std <- X_aug
  X_aug_std[, -1] <- scale(X_aug[, -1]) # Standardize (excluding Z)
  
  # Generate lambda sequence
  lambda_max <- max(abs(t(X_aug_std) %*% Y)) / n
  lambda_min <- lambda_max / 2e3
  k <- (0:(nlambda-1)) / nlambda
  lambda <- lambda_max * (lambda_min/lambda_max)^k
  
  # Fit gaussian model
  cv_fit <- glmnet::cv.glmnet(X_aug_std, Y, family = "gaussian", 
                      lambda = lambda, alpha = 1,
                      standardize = FALSE, standardize.response = FALSE)
  # Extract coefficients (remove intercept)
  coefs <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
  
  # Extract main effect coefficients
  Z_main <- coefs[2:(1 + p)]
  Z_main_tilde <- coefs[(2 + 2*p):(1 + 3*p)]
  
  # Extract interaction effect coefficients
  Z_int <- coefs[(2 + p):(1 + 2*p)]
  Z_int_tilde <- coefs[(2 + 3*p):(1 + 4*p)]
  
  # Compute W statistics with swap correction
  W_main <- (abs(Z_main) - abs(Z_main_tilde)) * (1 - 2*swap_main)
  W_int <- (abs(Z_int) - abs(Z_int_tilde)) * (1 - 2*swap_main)
  
  names(W_main) <- colnames(X)
  names(W_int) <- colnames(X)
  
  return(list(W_main = W_main, W_int = W_int))
}
