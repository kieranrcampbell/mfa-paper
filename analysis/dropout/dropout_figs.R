
library(tidyverse)
library(cowplot)

files <- dir("data", full.names = TRUE)
dfs <- lapply(files, read_csv)

df <- bind_rows(dfs)

df$pst_cor <- abs(df$pst_cor)
df$pst_cor_zi <- abs(df$pst_cor_zi)

names(df)[3:4] <- c("Normal", "Zero-inflated")

df <- gather(df, algorithm, correlation, -lambda, -rep, -pdrop)
df$lambda <- factor(as.character(df$lambda))

ggplot(df, aes(x = lambda, y = correlation, fill = algorithm)) + geom_boxplot() +
  xlab(expression(lambda)) + ylab("Correlation to true pseudotime") +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set1")

cor_plot <- last_plot()
  

ggplot(df, aes(x = lambda, y = 100 * pdrop)) + geom_boxplot() +
  xlab(expression(lambda)) + ylab("% zero expression")

cell_plot <- last_plot()

plot_grid(cor_plot, cell_plot, ncol = 1, labels = "AUTO")
