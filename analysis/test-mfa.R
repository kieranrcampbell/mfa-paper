library(ggplot2)
library(devtools)
library(tidyverse)
source("create_synthetic.R")

load_all("~/oxford/mfa/mfa")

synth <- create_synthetic(C = 100)

X <- synth$X; pst <- synth$pst; branch <- synth$branch


pca <- prcomp(t(X))

pca_tidy <- as_data_frame(pca$x[,1:2]) %>% 
  dplyr::mutate(pst, branch = factor(branch))

ggplot(pca_tidy, aes(x = PC1, y = PC2, colour = pst)) +
  geom_point() + scale_color_viridis()

ggplot(pca_tidy, aes(x = PC1, y = PC2, fill = branch)) +
  geom_point(shape = 21, alpha = 0.5) + 
  scale_fill_brewer(palette = "Set1")

y <- scale(t(X))
m <- mfa(y, iter = 2e4, b = 2, thin = 10, pc_initialise = 3, 
         alpha = 1e-2, beta = 1e-2,
         alpha_chi = 1e-2, beta_chi = 1e-2, prop_collapse = 1)
plot_mfa_trace(m)
ms <- summary(m)
pca_tidy$mfa_pseudotime <- ms$pseudotime
pca_tidy$mfa_branch <- ms$branch
pca_tidy$branch_certainty <- ms$branch_certainty

qplot(pca_tidy$pst, pca_tidy$mfa_pseudotime)

ggplot(pca_tidy, aes(x = PC1, y = PC2, colour = mfa_pseudotime)) +
  geom_point() + scale_color_viridis()

ggplot(pca_tidy, aes(x = PC1, y = PC2, fill = mfa_branch, size = branch_certainty)) +
  geom_point(shape = 21, alpha = 0.5) + 
  scale_fill_brewer(palette = "Set1")

ggplot(pca_tidy, aes(x = PC1, y = PC2, fill = branch_certainty)) +
  geom_point(shape = 21, alpha = 0.5) + 
  scale_fill_viridis()


# Test branch assignments -------------------------------------------------

# y <- scale(t(X))
cmap <- apply(m$traces$c_trace, 3, colMeans)
kmap <- apply(m$traces$k_trace, 3, colMeans)
tmap <- colMeans(m$traces$pst_trace)
taumap <- colMeans(m$traces$tau_trace)
chimap <- colMeans(m$traces$chi_trace)
plot(chimap)

# branch 1 probs
b1_probs <- sapply(seq_along(tmap), function(i) {
  dnorm(y[i, ], cmap[,1] + kmap[,1] * pst[i], 1 / sqrt(taumap), log = TRUE)
})
b2_probs <- sapply(seq_along(tmap), function(i) {
  dnorm(y[i, ], cmap[,2] + kmap[,2] * pst[i], 1 / sqrt(taumap), log = TRUE)
})
b1 <- colSums(b1_probs)
b2 <- colSums(b2_probs)

b12 <- b1 - b2

lse <- function(x) log(sum(exp(x - max(x)))) + max(x)
nr <- apply(cbind(b1, b2), 1, lse)
bp <- b1 - nr

qplot(nr, b1)
qplot(nr, b2)

pca_tidy$b12 <- b12
ggplot(pca_tidy, aes(x = PC1, y = PC2, fill = b12)) +
  geom_point(shape = 21) + 
  scale_fill_gradient2()

qplot(b1, b2, colour = pst, shape = factor(branch)) + 
  scale_color_viridis(name = "Pseudotime") +
  xlab("log p(X|branch1)") + ylab("log p(X|branch2)")


plot(b1_probs[,which.min(pst)])
plot(b2_probs[,which.min(pst)])


