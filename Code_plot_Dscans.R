# plot results from DSUITE Dinvestigate: various D-stat scans across genome

library(dplyr)
library(ggplot2)
library(ggExtra)


DATA="All_pops_Sbal_VCF.filter.99misd.maf5.SNPs"

setwd("/Users/johnburley/Dropbox (Brown)/Brown/Projects/Barnacle_speciation/")

P1="Iceland"
P2="Maine"
P3="RhodeIsland"
Outgroup="Canada"
SCAFFOLD="00030"
WINDOW=100
STEP=50

# could automate this readin with paste0:
data <- read.delim(paste0("Results/",DATA,"/Iceland_Maine_RhodeIsland_localFstats_Sbal_allSamples_allsites_100_50.txt"),
                    header=TRUE,
                    sep = "\t")

#data$D <- as.numeric(data$D)
#data$f_d <- as.numeric(data$f_d)
#data$f_dM <- as.numeric(data$f_dM)


# plot a scan of 
dat_scaf <- data %>%
  filter(chr == paste0("scaffold",SCAFFOLD)) %>%
  arrange(windowStart)

p <- ggplot(dat_scaf, aes(windowStart/1e6, f_dM)) +
  geom_point() +
  geom_smooth(method = "loess", span = 0.1) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = mean(dat_scaf$f_dM), linetype="dashed") +
  xlab("position along scaffold (Mbp)")
  
p1 <- ggMarginal(p, type="histogram")
p1

# explanation of plot
# y-axis histogram shows distribution of D stat
# x-axis histogram show the distribution of SNPs across the genome, seeing as window size is defined by 100 SNPs
# f_dM (Malinsky et al., 2015)

# is the distribution significantly differnt to zero?
t.test(dat_scaf$f_dM)

dev.copy2pdf(width = 7, height = 4, out.type = "pdf", 
             file = paste0("Results/",DATA,"/Figs/Iceland_Maine_RhodeIsland_localFstats_Sbal_allSamples_allsites_100_50.scaffold",SCAFFOLD,".pdf"))




