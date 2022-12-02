library(readr)
library(dplyr)
library(tidyr)
library(harmonicmeanp)
library(ggplot2)
library(scales)

Ns <- 1:15
is <- 0:999
path_1 <- 'subsampled_associations/'
path_2 <- '/assoc_ctl990_aff10'

read_length <- 150
tx_length <- 668
n_txs <- 4
n_aff <- 10
expected_overlaps_per_read <- (n_aff * (read_length - 1))/(n_txs * tx_length)

df <- data.frame(expected_overlap=Ns * expected_overlaps_per_read)

bonferroni_threshold <- 0.05 / 20000
y <- c()
for (N in Ns) {
  acc <- 0
  for (i in is) {
    path <- paste0(path_1, i, '/', N, path_2)
    df_ <- read_delim(path, col_names=c('kmer', '_', 'pval'), delim=" ")
    if (p.hmp(df_$pval) < bonferroni_threshold) {
      acc <- acc + 1
    }
  }
  y <- c(y, acc / length(is))
}

df$success_prop <- y
df <- rbind(data.frame(expected_overlap=0, success_prop=0), df)

breaks = c(0, 0.25, 0.5, 0.75, 1)

df %>% ggplot(aes(x = expected_overlap, y = success_prop)) +
        geom_line(size=1.5, colour='#00bda4') +
        # geom_smooth(se=F, aes(colour="Control: 990\nAffected: 10")) +
        scale_colour_manual(name="Cohort", values=c("#406661")) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0, 1)) +
        scale_x_continuous(limits=c(0,8)) +
        geom_hline(aes(yintercept=1), color="red", linetype="dashed") +
        labs(y = "Proportion of signals detected", x = "Expected # reads overlapping identifying region") +
        theme(legend.position = c(0.8, 0.2),
              text = element_text(family="Helvetica"),
              legend.title = element_text(size=20),
              plot.title = element_text(size = 30),
              axis.title = element_text(size=20),
              axis.text = element_text(size=12),
              legend.text = element_text(size=15))
