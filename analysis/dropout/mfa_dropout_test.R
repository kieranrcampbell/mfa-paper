
## Invoke via
## Rscript mfa_dropout_test.R [lambda] [rep]

library(devtools)
library(tidyverse)

load_all("~/mfa/mfa")
source("../create_synthetic.R")

args <- commandArgs(trailingOnly = TRUE)

lambda <- as.numeric(args[1])

rep <- as.numeric(args[2])

output_filename <- paste0("dropout_test_", lambda, "_", rep, ".csv")
output_file <- file.path("data", output_filename)

message(paste("Saving results to", output_file))

C <- 200
G <- 50

synth <- create_synthetic(C = C, G = G, lambda = lambda, model_dropout = TRUE, zero_negative = TRUE)

# X is gene-by-cell (like an ExpressionSet)
X <- synth$X; pst <- synth$pst

row_vars <- matrixStats::rowVars(X)
X <- X[row_vars > 0, ] # filter out zero variance genes
G <- nrow(X)

pc12 <- prcomp(t(X))$x[,1:2]
pc_cors <- apply(pc12, 2, cor, pst)
pc_initialise <- which.max(abs(pc_cors))

m <- mzi<- NULL

pst_cor <- pst_cor_zi <- NULL

if(any(X != 0)) {
	 m <- mfa(t(X), iter = 60000, thin = 30, 
             zero_inflation = FALSE,
             pc_initialise = pc_initialise,
             b = 2, tau_eta = 1e-5, tau_c = 0.1, tau_theta = 1e-2)

	     mzi <- mfa(t(X), iter = 60000, thin = 30, 
           zero_inflation = TRUE, 
           pc_initialise = pc_initialise,
           b = 2, tau_eta = 1e-5, tau_c = 0.1, tau_theta = 1e-2)

	   ms <- summary(m)
	   mzis <- summary(mzi)

	   pst_cor <- cor(pst, ms$pseudotime)
	   pst_cor_zi <- cor(pst, mzis$pseudotime)
} else {
  pst_cor <- rep(NA, ncol(X))
  pst_cor_zi <- rep(NA, ncol(X))
}

results <- data_frame(lambda, rep, pst_cor, pst_cor_zi, pdrop = mean(X == 0))

write_csv(results, output_file)