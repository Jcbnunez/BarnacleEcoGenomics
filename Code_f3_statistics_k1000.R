library(tidyverse)

a <- read_tsv("../reports/f3_statistics_k1000.txt")
b <- a %>% 
	mutate(p_val = pnorm(z_score)) %>% 
	mutate(p_val_adj = p.adjust(p_val)) %>% 
	mutate(log10_p_val = log10(p_val)) %>% 
	mutate(log10_p_val_adj = log10(p_val_adj))

write_tsv(b, "../reports/f3_statistics_k1000_adj.txt")

b <- b %>% mutate(index = row_number()) %>% 
	mutate(xlab = paste("", admix_pop, "; ", source_pop_1, ", ", source_pop_2, "", sep=""))

ggplot(b) + 
	geom_point(aes(index, f3, color=admix_pop)) + 
	geom_errorbar(aes(x=index, ymin=f3-f3_sd, ymax=f3+f3_sd, color=admix_pop)) + 
	labs(color="Target Population") + ylab("f3(Target; Source1, Source2)") + 
	scale_x_discrete(name="", limits=paste(b$xlab, sep="")) + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	geom_hline(yintercept=0, color="grey50", linetype="dashed") + 
	coord_flip() + theme(aspect.ratio=2)
ggsave("../reports/f3_statistics_k1000.pdf", scale=1.5)
