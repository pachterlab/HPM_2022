library(readr)
library(dplyr)
library(tidyr)
library(harmonicmeanp)
library(ggplot2)
library(scales)

Ns <- c(6, 10, 14, 18)
is <- 0:999
noise <- c('0.010000', '0.020000', '0.030000', '0.040000', '0.050000', '0.060000','0.070000','0.080000','0.090000','0.100000')
path_1 <- 'noise_subsampled_associations/'
path_2 <- '/assoc_ctl990_aff10'


read_length <- 150
tx_length <- 668
n_txs <- 4
n_aff <- 10
expected_overlaps_per_read <- (n_aff * (read_length - 1))/(n_txs * tx_length)

bonferroni_threshold <- 0.05 / 20000

df <- data.frame()
for (N in Ns) {
  for (n in noise) {
    acc <- 0
    for (i in is) {
      path <- paste0(path_1, i, '/', N, '/noise', n, path_2)
      df_ <- read_delim(path, col_names=c('kmer', '_', 'pval'), delim=" ")
      if (p.hmp(df_$pval) < bonferroni_threshold) {
        acc <- acc + 1
      }
    }
    df <- rbind(df, data.frame(expected_overlap=N * expected_overlaps_per_read, noise=as.numeric(n), success=acc/length(is)))

  }
}
df$expected_overlap <- as.factor(format(round(df$expected_overlap, 2), nsmall=2))

write_delim(df, 'success_proportion_noise.tsv', delim='\t')

# df <- read_delim('success_proportion_noise.tsv', delim="\t")

df %>% ggplot(aes(x = noise, y = success, colour = expected_overlap)) +
           geom_line(size=1.5) +
           scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0, 1)) +
           scale_x_continuous(limits=c(0.01,0.1)) +
           geom_hline(aes(yintercept=1), color="red", linetype="dashed") +
           labs(y = "Proportion of signals detected", x = "Magnitude of technical noise", colour="Expected # reads") +
           theme(legend.position = c(0.7, 0.3),
                 text = element_text(family="Helvetica"),
                 legend.title = element_text(size=20),
                 plot.title = element_text(size = 30),
                 axis.title = element_text(size=20),
                 axis.text = element_text(size=12),
                 legend.text = element_text(size=15))
