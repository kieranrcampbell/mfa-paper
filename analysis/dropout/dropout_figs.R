
library(tidyverse)

files <- dir("data", full.names = TRUE)
dfs <- lapply(files, read_csv)

df <- bind_rows(dfs)

df$pst_cor <- abs(df$pst_cor)
df$pst_cor_zi <- abs(df$pst_cor_zi)

names(df)[3:4] <- c("normal", "zero-inflated")

df <- gather(df, algorithm, correlation, -lambda, -rep)
df$lambda <- factor(as.character(df$lambda))

ggplot(df, aes(x = lambda, y = correlation, fill = algorithm)) + geom_boxplot()
