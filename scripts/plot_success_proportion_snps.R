library(readr)
library(dplyr)
library(tidyr)
library(harmonicmeanp)
library(ggplot2)
library(scales)
library(reshape)

path <- '/Users/kreldjarn/Dropbox/caltech/caltech_research/abundanceDBG_paper/success_proportion.tsv'

snps <- c(0.01, 0.05, 0.1, 0.5)

df <- read_delim(path, delim="\t")

molten <- melt(as.data.frame(df), id.vars="expected_overlap")

for (snp in snps) {
    molten[molten$variable==paste0("SNP rate ", snp),"expected_overlap"] <- molten[molten$variable==paste0("SNP rate ", snp),"expected_overlap"] / (1 - (snp + snp*snp + snp*snp*snp))
}

names(molten) <- c('expected_overlap', 'value', 'Genetic noise')

breaks = c(0, 0.25, 0.5, 0.75, 1)

molten %>% ggplot(aes(x = expected_overlap, y = value, colour = variable)) +
           geom_line(size=1.5) +
           scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0, 1)) +
           scale_x_continuous(limits=c(0,18)) +
           geom_hline(aes(yintercept=1), color="red", linetype="dashed") +
           labs(y = "Proportion of signals detected", x = "Expected # reads overlapping identifying region", colour="Genetic noise level") +
           theme(legend.position = c(0.7, 0.3),
                 text = element_text(family="Helvetica"),
                 legend.title = element_text(size=20),
                 plot.title = element_text(size = 30),
                 axis.title = element_text(size=20),
                 axis.text = element_text(size=12),
                 legend.text = element_text(size=15))
