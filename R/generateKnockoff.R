generateKnockoff <- function(X, Theta, Beta, seed = NULL) {
  D=nrow(Theta); V=ncol(Beta); N=rowSums(X)
  if(V == ncol(X)){
    colnames(Beta) <- colnames(X)
  }else if(V < ncol(X)){
    # Get the names of columns that are missing in Beta
    missing_cols <- setdiff(colnames(X), colnames(Beta))

    # For each missing column, add a column of zeros to Beta
    for (col in missing_cols) {
      Beta <- cbind(Beta, 0)
      colnames(Beta)[ncol(Beta)] <- col
    }

    # Reorder the columns of Beta to match the order in X
    Beta <- Beta[, colnames(X)]
  }
  # generate 1 sample
  generateSample <- function(N, theta, beta) {

    sample <- vector(length = N)
    z_d <- vector(length = N)
    for (n in 1:N) {
      z_n <- rmultinom(1, 1, theta)
      w_n <- rmultinom(1, 1, beta[which(z_n == 1),])
      sample[n] <- colnames(beta)[which(w_n == 1)]
      z_d[n] <- which(z_n == 1)
      names(z_d)[n] <- sample[n]
    }
    return(list(sample = sample, z = z_d))
  }

  # generate n samples
  cohort <- vector(mode = "list", length = D)
  z <- vector(mode = "list", length = D)
  for (d in 1:D) {
    sample.d <- generateSample(N[d], Theta[d, ], Beta)
    cohort[[d]] <- sample.d[["sample"]]
    z[[d]] <- sample.d[["z"]]
  }

  # collapse list of vectors into list of count tables
  sampleTaxaFreq <- lapply(cohort, table)

  # count matrix
  x_tilde <- matrix(data = 0, nrow = D, ncol = ncol(X))
  #rownames(x_tilde) <- rownames(Theta)
  colnames(x_tilde) <- colnames(Beta)

  for (d in 1:D) {
    x_tilde[d, names(sampleTaxaFreq[[d]])] <- sampleTaxaFreq[[d]]
  }

  return(x_tilde)
}
