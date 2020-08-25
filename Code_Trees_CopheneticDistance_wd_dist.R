#Packages
library(ape)
library(phylobase)
library(stringr)
library(magrittr)
library(ggpubr)
library(tidyverse)
library(forcats)
library(paleotree)
library(gmodels)

reverseList = function (list, simplify = FALSE) #from paleotree
{
    if (length(unique(sapply(list, length))) != 1) {
        stop("Not all lists equally long")
    }
    list1 <- list()
    for (i in 1:length(list[[1]])) {
        list1[[i]] <- sapply(list, function(x) x[[i]], simplify = simplify)
    }
    names(list1) <- names(list[[1]])
    return(list1)
}


#Load trees
filenames <- list.files("btrees/")
bootlist <- lapply(as.list(filenames), function(name){read.tree(paste("btrees/", name, sep=""))})
names(bootlist) <- str_replace_all(filenames, "\\.fasta\\.boottrees", "")
ngenes <- length(bootlist)
nboot <- length(bootlist[[1]])

#Homogenize labels
bootlist <- lapply(bootlist, function(item){lapply(item, function(t){t$tip.label <- str_replace_all(t$tip.label, "\\.\\w+\\.\\w",""); return(t)})})

#Make tree list + root + remove outgroup
bootlist <- lapply(bootlist, function(item){lapply(item, function(t){root.phylo(t, outgroup = "CAR")})})
bootlist <- lapply(bootlist, function(item){lapply(item, function(t){drop.tip(t, which(t$tip.label == "CAR"))})})

#branch length distances
Dists <- lapply(bootlist, function(item){lapply(item, cophenetic.phylo)}) 
Dists <- lapply(Dists, function(item){lapply(item, function(D){diag(D) <- rep(NA, length(diag(D))); return(D)})}) 
names(Dists) <- names(bootlist)


#Get BL distances within and between populations
pops <- unique(bootlist[[1]][[1]]$tip.label)
within <- lapply(pops, function(P){lapply(Dists, function(gene){lapply(gene, function(D){mean(D[which(rownames(D) == P),which(colnames(D) == P)], na.rm = T)})})})
names(within) <- pops
within <- lapply(reverseList(within), function(boot){reverseList(boot) %>% lapply(function(x){unlist(x) %>% mean()})})
between <- lapply(pops, function(P){lapply(Dists, function(gene){lapply(gene, function(D){mean(c(D[which(rownames(D) == P),which(colnames(D) != P)], D[which(rownames(D) != P),which(colnames(D) == P)]), na.rm = T)})})})
between <- lapply(reverseList(between), function(boot){reverseList(boot) %>% lapply(function(x){unlist(x) %>% mean()})})

#Determine balancing selection
bootBSnum <- lapply(as.list(1:ngenes), function(i){lapply(as.list(1:nboot), function(j){(within[[i]][[j]] - between[[i]][[j]])})})
names(bootBSnum) = names(bootlist)
bootBS <- lapply(bootBSnum, function(gene){lapply(gene, function(N){N>0})})

#plot balancing selection for bootstrap variants
bootBSnum <- data.frame(matrix(unlist(bootBSnum), ncol=length(bootBSnum)))
names(bootBSnum) <- names(bootlist)

#collect summary statistics
estimates = sapply(bootBSnum, function(x) t.test(x, mu=0, alternative = "greater")$estimate)  
stderr = sapply(bootBSnum, function(x) t.test(x, mu=0, alternative = "greater")$stderr)  
p_val = sapply(bootBSnum, function(x) t.test(x, mu=0, alternative = "greater")$p.value)  
ci_min = sapply(bootBSnum, function(x) t.test(x, mu=0, alternative = "greater")$conf.int[1])  
q25 = sapply(bootBSnum, function(x) quantile(x, 0.25))  
q75 = sapply(bootBSnum, function(x) quantile(x, 0.75))  

#make dataframe
data.frame(estimates,stderr,p_val, ci_min,q25,q75) %>% mutate(gene.etc = row.names(.)) %>% separate(gene.etc, into = c("geneid","etc"), sep = "\\.") %>% mutate(p_adjust = p.adjust(p_val, "bonferroni")) %>% mutate(signif = ifelse(.$p_adjust < 0.00000001, "sig","notsig") ) -> cophen_test

#save dataframe
write.table(cophen_test, file = "./cophen_test.txt", sep = "\t",quote = F ,row.names = F, col.names = T )

#plot graph
cophen_test %>% ggplot(aes(x=fct_reorder(geneid,estimates), y = estimates, ymin = q25, ymax=q75, fill = signif)) + geom_hline(yintercept = 0, size = 1 ) + geom_errorbar() + geom_point(size = 1.5, shape =21)  + theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_fill_manual(values = c("grey","steelblue3")) + theme(legend.pos = "none") -> cophne_ests

#save graph
ggsave("cophne_ests.pdf",cophne_ests, width = 3, height = 3)
