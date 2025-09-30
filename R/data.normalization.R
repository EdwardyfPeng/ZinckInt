###### Data normalization options ######
convert_to_proportions <- function(X) {
  if (!is.matrix(X)) {
    stop("Input must be a matrix.")
  }
  
  row_sums <- rowSums(X)
  
  if (any(row_sums == 0)) {
    warning("Some rows have zero total counts. These will result in NaN values.")
  }
  
  X_prop <- sweep(X, 1, row_sums, FUN="/")
  return(X_prop)
}

log_normalize <- function(X) {
  # Ensure input is a matrix
  if (!is.matrix(X)) {
    stop("Input must be a matrix.")
  }
  # Ensure all elements in the matrix are numeric
  if (!all(is.numeric(X))) {
    stop("Matrix contains non-numeric values.")
  }
  # Add a pseudo-count of 0.5 to each entry
  X <- X + 0.5
  # Make compositional by dividing each value by its row sum
  X <- sweep(X, 1, rowSums(X), FUN="/")
  # Apply log transformation
  Z <- log(X)
  return(Z)
}

CLR_normalize <- function(X){
   # Ensure input is a matrix
  if (!is.matrix(X)) {
    stop("Input must be a matrix.")
  }
  # Ensure all elements in the matrix are numeric
  if (!all(is.numeric(X))) {
    stop("Matrix contains non-numeric values.")
  }
  # Add a pseudo-count of 0.5 to each entry
  X <- X + 0.5
  # Make compositional by dividing each value by its row sum
  X <- sweep(X, 1, rowSums(X), FUN="/")
  # apply CLR (centered log-ratio transformation) to proportions
  mu   <- rowMeans(log(X))                        
  clrX <- sweep(log(X), 1, mu, FUN = "-") 
  return(clrX)
}
