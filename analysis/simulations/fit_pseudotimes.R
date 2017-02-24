
#'
#' Invoke via
#' Rscript fit_pseudotimes.R [n_cells] [n_genes] [p_transient] [rep]
#' and output results will be written to
#' data/pseudotime_results_[n_cells]_[n_genes]_[p_transient]_[rep].csv

library(dpt)
library(monocle)
library(devtools)
library(tidyverse)

load_all("oxford/mfa/mfa")
source("../create_synthetic.R")

fit_dpt_pseudotime <- function(X, pst) {
  root_cell <- which.min(pst)
  ts <- Transitions(t(X))

  pt <- tryCatch(dpt(ts, branching = TRUE, root = root_cell)$DPT,
                 error = function(e) {
                    message("DPT failed")
                   return(rep(NA, ncol(X)))
                 })
  return(pt)
}

fit_monocle_pseudotime <- function(X, pst) {
  root_cell <- which.min(pst)
  cds <- newCellDataSet(X)
  cds <- reduceDimension(cds, norm_method = "none")
  cds <- orderCells(cds, num_paths = 2)
  cds <- orderCells(cds, num_paths = 2, root_state = cds$State[root_cell])
  return(cds$Pseudotime)
}

fit_tscan_pseudotime <- function(X, pst) {
  eclust <- exprmclust(X)
  TSCANorder(eclust, orderonly = FALSE)
}

fit_mfa_pseudotime <- function(X, pst) {
  pc12 <- prcomp(t(X))$x[,1:2]
  pc_cors <- apply(pc12, 2, cor, pst)
  pc_initialise <- which.max(abs(pc_cors))

  m <- mfa(t(X), iter = 6000, b = 2, pc_initialise = pc_initialise,
           alpha = 0.1, beta = 0.1)
  ms <- summary(m)
  return(ms$pseudotime)
}


args <- commandArgs(trailingOnly = TRUE)

# print(args)

C <- as.numeric(args[1]) # Number of genes
G <- as.numeric(args[2]) # Number of genes
p_transient <- as.numeric(args[3]) # Proporiton of transient genes
rep <- as.numeric(args[4])

output_name <- paste("pseudotime", "results", C, G, p_transient, rep, sep = "_")
output_file <- file.path("data", paste0(output_name, ".csv"))

synth <- create_synthetic(C = C, G = G, p_transient = p_transient)

# X is gene-by-cell (like an ExpressionSet)
X <- synth$X; pst <- synth$pst

output_df <- data_frame(
  monocle = fit_monocle_pseudotime(X, pst),
  dpt = fit_dpt_pseudotime(X, pst),
  mfa = fit_mfa_pseudotime(X, pst),
  true_pst = pst
)

output_df <- mutate(output_df,
                    C = C, G = G, p_transient = p_transient, rep = rep)

write_csv(output_df, output_file)
