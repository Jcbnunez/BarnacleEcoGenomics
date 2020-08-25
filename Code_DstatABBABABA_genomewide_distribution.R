# investigate distribution of D (ABBABABA) stats across scaffolds.

# possible extension could be to merge this with genome annotation. compare coding vs non coding etc. 

library(dplyr)
library(ggplot2)
library(ggExtra)

DATA="All_pops_Sbal_VCF.filter.99misd.maf5.SNPs"

setwd("/Users/johnburley/Dropbox (Brown)/Brown/Projects/Barnacle_speciation/")

P1="Iceland"
P2="Maine"
P3="RhodeIsland"
Outgroup="Canada"
WINDOW=100
STEP=50

# could automate this readin with paste0:
data <- read.delim(paste0("Results/",DATA,"/Iceland_Maine_RhodeIsland_localFstats_Sbal_allSamples_allsites_100_50.txt"),
                   header=TRUE,
                   sep = "\t")

data_scaf_summary <- data %>%
  group_by(chr) %>%
  dplyr::summarise(n.windows = n(),
                   f_dM.mean = mean(f_dM)) %>%
    filter(n.windows > 6)

hist(data_scaf_summary$f_dM.mean)

summary(data_scaf_summary$f_dM.mean)

plot(data_scaf_summary$n.windows, data_scaf_summary$f_dM.mean, label = data_scaf_summary$chr)

ggplot(data_scaf_summary, aes(n.windows,f_dM.mean)) +
  #geom_point() +
  geom_text(aes(label=gsub("scaffold","", chr)))

