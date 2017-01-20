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

C <- 100 # cells
G <- 40 # genes

p_transient <- 0.5 # proportion of gene behaviour that should be transient

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
X[transient_genes_common,] <- t(sapply(transient_genes_common, function(g) {
  scale <- rlnorm(1, log(0.05), 0.5)
  reverse <- sample(c(T,F), 1)
  mu <- 2 * phi[g, 1] * transient(pst, scale = scale, reverse = reverse)
  rnorm(length(mu), mu, gsd[g])
}))

# Deal with bifurcating ones
X[transient_genes_bifurcating,] <- t(sapply(transient_genes_bifurcating, function(g) {
  which_nonzero <- which(k[g,] != 0) # we're going to make this one transient
  scale <- rlnorm(1, log(0.05), 0.3)
  reverse <- k[g, which_nonzero] < 0
  mu <- phi[g, which_nonzero] * transient(pst, location = 0.75, 
                                              scale = scale, reverse = reverse)
  
  cells_on_constant_branch <- which(branch != which_nonzero)
  cells_on_transient_branch <- which(branch == which_nonzero)
  
  y <- rep(NA, C)
  y[cells_on_transient_branch] <- rnorm(length(cells_on_transient_branch), 
                                        mu[cells_on_transient_branch], gsd[g])
  y[cells_on_constant_branch] <- X[g, cells_on_constant_branch]
  return( y )
}))

stop("done")


# Visualise expression ----------------------------------------------------

x_tidy <- t(X) %>% as_data_frame() %>% 
  mutate(pst, branch = factor(branch)) %>% 
  gather(gene, expression, -pst, -branch)

ggplot(x_tidy, aes(x = pst, y = expression, fill = branch)) + 
  geom_point(shape = 21, alpha = 0.5) + 
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_brewer(palette = "Set1")



# 
# fname <- file.path("data", "synthetic.h5")
# if(!file.exists(fname)) {
#   h5createFile(fname)
# }
# 
# h5createGroup(fname, "basic_branching") # we'll call this dataset basic_branching
# h5write(t(X), fname, "basic_branching/X")
# h5write(pst, fname, "basic_branching/pseudotime")
# h5write(branch, fname, "basic_branching/branch_assignment")
# h5write(k, fname, "basic_branching/k")
# h5write(phi, fname, "basic_branching/phi")
# h5write(delta, fname, "basic_branching/delta")
# h5write(gsd, fname, "basic_branching/stdev")
# 

pca <- prcomp(t(X))

pca_tidy <- as_data_frame(pca$x[,1:2]) %>% 
  mutate(pst, branch = factor(branch))

ggplot(pca_tidy, aes(x = PC1, y = PC2, colour = pst)) +
  geom_point() + scale_color_viridis()

ggplot(pca_tidy, aes(x = PC1, y = PC2, fill = branch)) +
  geom_point(shape = 21, alpha = 0.5) + 
  scale_fill_brewer(palette = "Set1")

pca_tidy$mfa_pst <- ms$pseudotime
pca_tidy$mfa_branch <- ms$branch
pca_tidy$branch_certainty <- ms$branch_certainty

ggplot(pca_tidy, aes(x = PC1, y = PC2, colour = mfa_pst)) +
  geom_point() + scale_color_viridis()

ggplot(pca_tidy, aes(x = PC1, y = PC2, colour = mfa_branch, size = 1 / log(branch_certainty))) +
  geom_point() 


# 
# ggsave("figs/basic_branching.png", width = 10, height = 8)
