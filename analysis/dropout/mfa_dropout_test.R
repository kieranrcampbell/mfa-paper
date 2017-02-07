
## Invoke via
## Rscript mfa_dropout_test.R [lambda] [rep]

library(devtools)
library(tidyverse)

load_all("~/mfa/mfa")
source("../create_synthetic.R")

args <- commandArgs(trailingOnly = TRUE)

lambda <- as.numeric(args[1])
rep <- as.numeric(args[2])

output_filename <- paste0("dropout_test_", lambda, "_", rep, ".Rdata")
output_file <- file.path("data", output_filename)

C <- 200
G <- 50

synth <- create_synthetic(C = C, G = G, lambda = lambda, model_dropout = TRUE)

# X is gene-by-cell (like an ExpressionSet)
X <- synth$X; pst <- synth$pst

pc12 <- prcomp(t(X))$x[,1:2]
pc_cors <- apply(pc12, 2, cor, pst)
pc_initialise <- which.max(abs(pc_cors))

m <- NULL


m <- mfa(t(X), iter = 20000, thin = 10, 
             zero_inflation = TRUE, lambda = lambda,
             pc_initialise = pc_initialise,
             b = 2, tau_eta = 1e-5, tau_c = 0.1, tau_theta = 1e-2)

mzi <- mfa(t(X), iter = 20000, thin = 10, 
           zero_inflation = TRUE, lambda = lambda,
           pc_initialise = pc_initialise,
           b = 2, tau_eta = 1e-5, tau_c = 0.1, tau_theta = 1e-2)


save(m, mzi, synth, lambda, model_dropout, rep, file = output_file)