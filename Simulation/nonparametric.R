library(rstan)
library(glmnet)
library(knockoff)    
library(topicmodels) 
library(randomForest)
library(tidyverse)
library(ggplot2)
library(stats)
library(ranger)
library(kosel)
library(reshape2)
library(lubridate)
library(Matrix)

load("count.Rdata")

####### Stan code for fitting Zinck model #######
zinck_code <- "data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num words
  int<lower=0> D; // num docs
  int<lower=0> n[D, V]; // word counts for each doc

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma1;
  vector<lower=0>[V] gamma2;
  vector<lower=0, upper=1>[V] delta;
}
parameters {
  simplex[K] theta[D]; // topic mixtures
  vector<lower=0,upper=1>[V] zeta[K]; // zero-inflated betas
}


transformed parameters {
  vector<lower=0>[V] beta[K];
  for (k in 1:K) {
	beta[k,1] =  zeta[k,1];
  for (m in 2:V) {
    beta[k,m] = zeta[k,m]*prod(1 - zeta[k,1:(m - 1)]);  // stick breaking
  }
  }
  for (k in 1:K) {
      beta[k]=beta[k]/sum(beta[k,1:V]);  // GD construction
  }
}


model {
  for (d in 1:D) {
    theta[d] ~ dirichlet(alpha);  
  }

  for (k in 1:K) {
    for (m in 1:V) {
      if (zeta[k,m]==0){  // Zero-inflated beta likelihood
        target += bernoulli_lpmf(1 | delta[m]);
      }else{
        target += bernoulli_lpmf(0 | delta[m]) + beta_lpdf(zeta[k,m] | gamma1[m], gamma2[m]);
      }
		}
  }

  for (d in 1:D) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];
    for (k in 2:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    eta = eta/sum(eta[1:V]);
    n[d] ~ multinomial(eta);  // generation of each sample
  }
}
"
stan.model = stan_model(model_code = zinck_code)

X_sim <- generate_data_int(p=200, seed=1)$X
Y_sim <- generate_data_int(p=200, seed=1)$Y
Z_sim <- generate_data_int(p=200, seed=1)$Z
true_main <- generate_data_int(p=200, seed = 1)$S_beta
true_int <- generate_data_int(p=200, seed=1)$S_gamma

### ======== Apply differnt methods to do feature selection(set target_fdr=0.2, offset=1(knockoff+)) ======== ###

## Method1: Zinck ##
dlt <- c()
for(t in (1:ncol(X_sim))){
   dlt[t] <- 1-mean(X_sim[,t]>0)
}

zinck_stan_data <- list(
  K = 15,
  V = ncol(X_sim),
  D = nrow(X_sim),
  n = X_sim, 
  alpha = rep(0.2, 15),  # change from 1.0 to 0.2 for more zero-inflated knockoff (0.1 will genearte error)
  gamma1 = rep(0.5, ncol(X_sim)),
  gamma2 = rep(10.0, ncol(X_sim)),
  delta = dlt
)

fit1 <- vb(stan.model,
           data = zinck_stan_data,
           algorithm = "meanfield",
           importance_resampling = TRUE,
           iter = 10000,
           tol_rel_obj = 0.01, 
           elbo_samples = 500
)

theta <- fit1@sim[["est"]][["theta"]]
beta <- fit1@sim[["est"]][["beta"]]

Xsim_tilde <- generateKnockoff(X_sim, theta, beta, seed=1) ## getting the knockoff copy

index_est_zinck_her <- suppressWarnings(zinck.filter.interaction(X = X_sim,
                                                                 X_tilde = Xsim_tilde,
                                                                 Y = Y_sim,
                                                                 Z = Z_sim,
                                                                 model = "glmnet",
                                                                 fdr = 0.2,
															     offset = 1,
															     seed = 1,
                                                                 heredity = "strong"
                                                                 ))
index_est_zinck_sep <- suppressWarnings(zinck.filter.interaction(X = X_sim,
                                                                 X_tilde = Xsim_tilde,
                                                                 Y = Y_sim,
                                                                 Z = Z_sim, 
																 model = "glmnet"
                                                                 fdr = 0.2,
																 offset = 1,
																 seed = 1,
                                                                 heredity = "separate"
                                                                 ))

## Method2: MX-KF ##
Xlog <- log_normalize(X_sim)
Xlog_tilde <- create.second_order(Xlog)
index_est_mxkf_sep <- zinck.filter.interaction(X = Xlog, 
                                               X_tilde = Xlog_tilde, 
                                               Y = Y_sim, 
                                               Z = Z_sim,
                                               model="glmnet", 
                                               fdr=0.2, 
											   offset = 1,
											   seed = 1,
                                               heredity = "separate")
index_est_mxkf_her <- zinck.filter.interaction(X = Xlog, 
                                               X_tilde = Xlog_tilde, 
                                               Y = Y_sim, 
                                               Z = Z_sim,
                                               model = "glmnet", 
                                               fdr=0.2, 
											   offset = 1,
											   seed=1,
                                               heredity = "strong")

## Method3: LDA-KF ##
df.LDA <- as(as.matrix(X1),"dgCMatrix")
vanilla.LDA <- LDA(df.LDA,k=8,method="VEM") 
theta.LDA <- vanilla.LDA@gamma
beta.LDA <- vanilla.LDA@beta
beta.LDA <- t(apply(beta.LDA,1,function(row) row/sum(row)))
Xsim_tilde.LDA <- generateKnockoff(X_sim,theta.LDA,beta.LDA,seed=1) ## Generating vanilla LDA knockoff copy
index_est_lda_h <- zinck.filter.interaction(X = X_sim,
                                            X_tilde = Xsim_tilde.LDA,
                                            Y = Y_sim,
                                            Z = Z_sim,
                                            model="glmnet",
                                            fdr=0.2,
											offset=1,
											seed=1, 
											heredity = "strong")
index_est_lda_sep <- zinck.filter.interaction(X = X_sim,
                                              X_tilde = X1_tilde.LDA,
                                              Y = Y_sim,
                                              Z = Z_sim,
                                              model="glmnet",
                                              fdr=0.2,
											  offset=1,
											  seed=1, 
											  heredity = "separate")

## Method4: DeepLINK (Python Code) ##
import DeepLINK as dl
from PCp1_numFactors import PCp1 as PCp1
import numpy as np
import pandas as pd
import keras
from keras.layers import Dense, Dropout
from keras.models import Sequential
from pairwise_connected_layer import PairwiseConnected
from itertools import combinations
from keras.callbacks import EarlyStopping
import tensorflow as tf
import random

count =  pd.read_csv("count.csv")
count = count.values

def generate_data_int(p, seed):
    np.random.seed(seed)
    # Order columns by decreasing abundance
    col_sums = np.sum(count, axis=0)
    sorted_indices = np.argsort(-col_sums)  # Descending order
    dcount = count[:, sorted_indices]

    # Filter features with >20% non-zero samples
    norm_count = count / np.sum(count, axis=1, keepdims=True)
    col_means = np.mean(norm_count > 0, axis=0)
    valid_indices = np.where(col_means > 0.2)[0]
    sorted_valid_indices = valid_indices[np.argsort(-col_means[valid_indices])]

    # Select top p features (p could be 200, 300, 400)
    dcount = count[:, sorted_valid_indices][:, :p]

    # Randomly sample 500 observations
    n_total = dcount.shape[0]
    sel_index = np.sort(np.random.choice(n_total, size=500, replace=False))
    dcount = dcount[sel_index, :]

    # Add pseudocount and calculate proportions
    original_OTU = dcount + 0.5
    seq_depths = np.sum(original_OTU, axis=1)
    Pi = original_OTU / seq_depths[:, np.newaxis]
    n = Pi.shape[0]

    Z = np.random.permutation(np.concatenate([np.ones(n // 2, dtype=int),
                                              -np.ones(n // 2, dtype=int)]))
    biomarker_main_idx = np.random.choice(np.arange(200), size=40, replace=False)
    biomarker_int_idx = np.random.choice(biomarker_main_idx, size=20, replace=False)

    # Average true proportion per feature
    Cj = Pi.mean(axis=0)

    # ----- Main effects (beta) -----
    beta = np.zeros(p)
    U_main = np.random.uniform(6, 12, size=40)
    sign_beta = np.random.choice([1, -1], size=40)
    beta_vals = sign_beta * (U_main / np.sqrt(Cj[biomarker_main_idx]))
    beta[biomarker_main_idx] = beta_vals

    # ----- Interaction effects (gamma) -----
    gamma = np.zeros(p)
    sign_gamma = np.random.choice([1, -1], size=20)
    gamma_vals = sign_gamma * (np.abs(beta[biomarker_int_idx]) / 2.0)
    gamma[biomarker_int_idx] = gamma_vals

    # ----- Simulate continuous outcome Y ~ N(mu, 1) -----
    alpha = 1.0
    fPi = Pi                                # shape (n, p) # remove quadratic terms
    mu_main = fPi @ beta                    # shape (n,)
    mu_int = (fPi @ gamma) * Z              # shape (n,)
    mu = alpha * Z + mu_main + mu_int
    Y = np.random.normal(loc=mu, scale=1.0, size=n)

    # ----- Generate final count data via multinomial sampling -----

    sim_count = np.zeros((n, p), dtype=int)
    for i in range(n):
        probs = Pi[i]
        sim_count[i] = np.random.multinomial(seq_depths[i], probs)
    
    return {
        "X": sim_count,
        "Y": Y,
        "Z": Z,
        "S_beta": biomarker_main_idx,
        "S_gamma": biomarker_int_idx,
        "beta": beta,
        "gamma": gamma
    }

aut_epoch = 100 # number of autoencoder training epochs
aut_loss = 'mean_squared_error' # loss function used in autoencoder training
aut_verb = 0 # verbose level of autoencoder
mlp_epoch = 100 # number of mlp training epochs
mlp_loss = 'binary_crossentropy'
#mlp_loss = 'mean_squared_error' # loss function used in mlp training
dnn_loss = 'binary_crossentropy'
dnn_verb = 0
aut_met = 'relu'
dnn_met = 'elu'
mlp_verb = 0 # verbose level of mlp
l1 = 0.001 # l1 regularization factor in mlp
lr = 0.001 # learning rate for mlp training

def run_deeplink_method(X, Y, Z, q):
    try:
        X1 = X.astype(np.float64)
        X1 -= np.mean(X1, axis=0)
        std_dev = np.std(X1, axis=0, ddof=1)
        std_dev[std_dev == 0] = 1
        X1 /= std_dev
        norms = np.sqrt(np.sum(X1 ** 2, axis=0))
        norms[norms == 0] = 1
        X2 = X1 / norms
        
        r_hat = PCp1(X2, 15)
        
        X_knockoff_combined = dl.knockoff_construct(X1, r_hat, 'elu', aut_epoch, aut_loss, aut_verb)
        p = X1.shape[1]
        
        X_original = X_knockoff_combined[:, :p]
        X_tilde = X_knockoff_combined[:, p:]
        
        X_Z = X_original * Z.reshape(-1, 1)
        X_tilde_Z = X_tilde * Z.reshape(-1, 1)
        
        es = EarlyStopping(monitor='val_loss', patience=30, verbose=0)
        
        Xnew_main = np.concatenate([X_original, X_tilde], axis=1)
        
        dp_main = Sequential()
        dp_main.add(PairwiseConnected(input_shape=(2 * p,)))
        dp_main.add(Dense(p, activation='elu', kernel_regularizer=keras.regularizers.l1(l1=l1)))
        dp_main.add(Dense(1, activation=None))
        dp_main.compile(loss=mlp_loss, optimizer=keras.optimizers.Adam(learning_rate=lr))
        dp_main.fit(Xnew_main, Y, epochs=mlp_epoch, batch_size=32, verbose=mlp_verb, 
                    validation_split=0.1, callbacks=[es])
        
        weights_main = dp_main.get_weights()
        w_main = weights_main[1] @ weights_main[3]
        w_main = w_main.reshape(p, )
        z_main = weights_main[0][:p]
        z_tilde_main = weights_main[0][p:]
        W_main = (w_main * z_main) ** 2 - (w_main * z_tilde_main) ** 2
        
        Xnew_int = np.concatenate([X_Z, X_tilde_Z], axis=1)
        
        dp_int = Sequential()
        dp_int.add(PairwiseConnected(input_shape=(2 * p,)))
        dp_int.add(Dense(p, activation='elu', kernel_regularizer=keras.regularizers.l1(l1=l1)))
        dp_int.add(Dense(1, activation=None))
        dp_int.compile(loss=mlp_loss, optimizer=keras.optimizers.Adam(learning_rate=lr))
        dp_int.fit(Xnew_int, Y, epochs=mlp_epoch, batch_size=32, verbose=mlp_verb, 
                   validation_split=0.1, callbacks=[es])
        
        weights_int = dp_int.get_weights()
        w_int = weights_int[1] @ weights_int[3]
        w_int = w_int.reshape(p, )
        z_int = weights_int[0][:p]
        z_tilde_int = weights_int[0][p:]
        W_int = (w_int * z_int) ** 2 - (w_int * z_tilde_int) ** 2
        
        # Separate
        selected_main_sep = dl.knockoff_select(W_main, q, ko_plus=False)
        selected_int_sep = dl.knockoff_select(W_int, q, ko_plus=False)
        
        # Hierarchical
        selected_main_her = dl.knockoff_select(W_main, q, ko_plus=False)
        if len(selected_main_her) > 0:
            W_int_restricted = W_int[selected_main_her]
            selected_indices = dl.knockoff_select(W_int_restricted, q, ko_plus=False)
            selected_int_her = selected_main_her[selected_indices] if len(selected_indices) > 0 else []
        else:
            selected_int_her = []
        
        return {
            'selected_main_sep': selected_main_sep,
            'selected_int_sep': selected_int_sep,
            'selected_main_her': selected_main_her,
            'selected_int_her': selected_int_her,
            'success': True
        }
    
    except Exception as e:
        print(f"Error in run_deeplink_method: {e}")
        return {
            'selected_main_sep': [],
            'selected_int_sep': [],
            'selected_main_her': [],
            'selected_int_her': [],
            'success': False
        }
