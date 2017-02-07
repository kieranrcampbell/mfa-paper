library(ggplot2)
library(viridis)

# set.seed(123L)

#' Sigmoid function for activations
sigmoid <- function(t, phi, k, delta) {
  return( 2 * phi / (1 + exp(-k*(t - delta))))
}

transient <- function(t, location = 0.5, scale = 0.01, reverse = FALSE) {
  y <- exp(- 1 / (2 * scale) * (t - location)^2)
  if(reverse) y <- 1 - y
  return(y)
}

create_synthetic <- function(C = 100, G = 40, p_transient = 0,
                             zero_negative = FALSE, model_dropout = FALSE,
                             lambda = 1) {
  
  # C <- 100 # cells
  # G <- 40 # genes
  
  # p_transient <- 1 # proportion of gene behaviour that should be transient
  
  branch <- rbinom(C, 1, 0.5)
  
  gsd <- sqrt(1 / rgamma(G, 2, 2))
  
  
  ## We assume first G / 2 (= 20) genes are common to both branches, and the 
  ## final G / 2 genes exhibit branching structure. We want to build in the 
  ## fact that delta < 0.5 for the common genes and delta > 0.5 for the 
  ## branch specific genes
  
  k <- replicate(2, runif(G, 5, 10) * sample(c(-1, 1), G, replace = TRUE))
  phi <- replicate(2, runif(G, 5, 10))
  delta <- replicate(2, runif(G, 0.5, 1))
  
  # Non bifurcating genes
  inds <- 1:(G / 2)
  
  # Bifurcating genes
  inds2 <- (G/2 + 1):G
  
  # For non-bifurcating genes, set behaviour identical across the two branches
  k[, 2] <- k[inds, 1]
  
  k[inds2, ] <- t(apply(k[inds2, ], 1, function(r) r * sample(c(0, 1))))
  
  phi[, 1] <- phi[, 2]
  delta[inds, 2] <- delta[inds, 1] <- runif(G / 2, 0, 0.5)
  
  ## Now make it look like a branching process
  for(r in inds2) {
    whichzero <- which(k[r,] == 0)
    nonzero <- which(k[r,] != 0)
    k_sign <- sign(k[r,nonzero])
    if(k_sign == 1) {
      phi[r, whichzero] <- 0
    } else {
      phi[r, whichzero] <- 2 * phi[r, nonzero]
    }
  }
  
  pst <- runif(C)
  
  X <- sapply(seq_along(branch), function(i) {
    k_i <- k[, branch[i] + 1]
    phi_i <- phi[, branch[i] + 1]
    delta_i <- delta[, branch[i] + 1]
    mu <- sigmoid(pst[i], phi_i, k_i, delta_i)
    rnorm(length(mu), mu, gsd)
  })
  
  ## Now let's add in the transient genes
  
  transient_genes <- sample(C, round(p_transient * C))
  transient_genes_common <- intersect(transient_genes, inds)
  transient_genes_bifurcating <- intersect(transient_genes, inds2)
  
  
  # Deal with non-bifurcating ones
  if(length(transient_genes_common) > 0) {
    X[transient_genes_common,] <- t(sapply(transient_genes_common, function(g) {
      scale <- rlnorm(1, log(0.05), 0.5)
      reverse <- sample(c(T,F), 1)
      mu <- 2 * phi[g, 1] * transient(pst, scale = scale, reverse = reverse)
      rnorm(length(mu), mu, gsd[g])
    }))
  }
  
  # Deal with bifurcating ones
  if(length(transient_genes_bifurcating) > 0) {
    X[transient_genes_bifurcating,] <- t(sapply(transient_genes_bifurcating, function(g) {
      which_nonzero <- which(k[g,] != 0) # we're going to make this one transient
      scale <- rlnorm(1, log(0.05), 0.3)
      reverse <- k[g, which_nonzero] < 0
      mu <- 2 * phi[g, which_nonzero] * transient(pst, location = 0.75, 
                                                  scale = scale, reverse = reverse)
      
      cells_on_constant_branch <- which(branch != which_nonzero)
      cells_on_transient_branch <- which(branch == which_nonzero)
      
      y <- rep(NA, C)
      y[cells_on_transient_branch] <- rnorm(length(cells_on_transient_branch), 
                                            mu[cells_on_transient_branch], gsd[g])
      y[cells_on_constant_branch] <- X[g, cells_on_constant_branch]
      return( y )
    }))
  }
  
  if(zero_negative) {
    X[X < 0] <- 0
  }
  
  if(model_dropout && lambda < Inf) {
    drop_probs <- t(apply(X, 1, function(x) exp(-lambda * x)))
    for(g in seq_len(G)) {
      drop <- runif(C) < drop_probs[g, ]
      X[g,drop] <- 0
    }
  }
  
  # row_vars <- matrixStats::rowVars(X)
  # if(any(row_vars == 0)) {
  #   zero_genes <- which(row_vars == 0)
  #   for(g in zero_genes) {
  #     N_cells <- sample(seq_len(C / 5), 1) # how many cells to add stochastic expression to
  #     which_cells <- sample(seq_len(C), N_cells)
  #     X[g,which_cells] <- runif(N_cells, 0, 0.5) # add small noise
  #   }
  # }
  
  list(X = X, branch = branch, pst = pst, k = k, phi = phi,
       delta = delta, p_transient = p_transient)
}
  
