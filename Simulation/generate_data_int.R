###### Data Generating Process ######
generate_data_int <- function(p, seed){
  # Ordering the columns with decreasing abundance
  dcount <- count[, order(decreasing = TRUE, colSums(count, na.rm = TRUE), 
                          apply(count, 2L, paste, collapse = ''))] 
  
  ## Randomly sampling 500 samples from 574 observations
  set.seed(seed)
  norm_count <- count / rowSums(count)
  col_means <- colMeans(norm_count > 0)
  indices <- which(col_means > 0.2)
  sorted_indices <- indices[order(col_means[indices], decreasing = TRUE)]
  
  dcount <- count[, sorted_indices][, 1:p] #top p (p = 200,300,or 400) most abundant features.
  sel_index <- sort(sample(1:nrow(dcount), 500))
  dcount <- dcount[sel_index, ]
  original_OTU <- dcount + 0.5
  seq_depths <- rowSums(original_OTU)
  Pi <- sweep(original_OTU, 1, seq_depths, "/")
  n <- nrow(Pi)
  
  ## Generate binary exposure variable Z with equal levels (+1/-1)
  Z <- sample(rep(c(1,-1), each = n/2))
  
  ## Calculate size effect
  # randomly select 40 biomarkers with nonzero main effects from the top 200 most abundant feature
  biomarker_main_idx <- sample(1:200, 40, replace = FALSE)
  
  # randomly select 20 biomarkers to have nonzero interaction effects
  biomarker_int_idx <- sample(biomarker_main_idx, 20, replace = FALSE)
  
  # Cj is the average true proportion of feature j
  Cj <- colMeans(Pi)
  
  # Main effects
  beta <- rep(0,p)
  U_main <- runif(40, min = 50, max = 100)
  sign_beta <- sample(c(1,-1), size = 40, replace = TRUE)
  beta_vals <- sign_beta * (U_main / sqrt(Cj[biomarker_main_idx]))
  beta[biomarker_main_idx] <- beta_vals
  
  # Interaction effects
  gamma <- rep(0, p)
  sign_gamma <- sample(c(1, -1), size = 20, replace = TRUE)
  gamma_vals <- sign_gamma * (abs(beta[biomarker_int_idx])/2)
  gamma[biomarker_int_idx] <- gamma_vals
  
  ## simulated the continuous outcome Yi from a normal distribution with unit variance.
  alpha <- 1
  fPi <- 0.5 * (Pi ^ 2) + Pi                 # n x p
  mu_main <- as.vector(fPi %*% beta)             # n
  mu_int <- as.vector(fPi %*% gamma) * Z        # n
  mu <- alpha * Z + mu_main + mu_int
  Y <- rnorm(n, mean = mu, sd = 1)
  
  ## Generate final count data using multinomial sampling
  
  # Sample sequencing depths from template data
  template_seq_depths <- rowSums(count)
  drawn_depths <- sample(template_seq_depths, size = n, replace = TRUE)
  
  # Generate simulated count data
  sim_count <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n) {
    sim_count[i, ] <- rmultinom(1, size = drawn_depths[i], prob = Pi[i, ])
  }
  colnames(sim_count) <- colnames(Pi) 
  
  return(list(
    X = sim_count,
    Y = Y,
    Z = Z,
    S_beta = biomarker_main_idx,
    S_gamma = biomarker_int_idx,
    beta = beta,
    gamma = gamma
  ))
}                     
