library(readr)
library(dplyr)
library(tidyr)
library(harmonicmeanp)

path_1 <- 'subsampled_associations_representative/snps0.01/'
path_2 <- '/25/assoc_ctl990_aff10'
is <- 0:999

bonferroni_threshold <- 0.05 / 20000


acc <- 0
for (i in is) {
  path <- paste0(path_1, i, path_2)
  df_ <- read_delim(path, col_names=c('kmer', '_', 'pval'), delim=" ")
  if (p.hmp(df_$pval) < bonferroni_threshold) {
    acc <- acc + 1
  }
}
print(acc/length(is))
