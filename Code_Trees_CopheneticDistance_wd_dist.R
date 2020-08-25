#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Packages
library(ape)
library(phylobase)
library(stringr)
library(magrittr)

reverseList = function (list, simplify = FALSE) 
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


#Tree files
tree_in <- read.tree(args[1])
gene=as.character(args[2])

#Homogenize labels
tree_in$tip.label <- str_replace_all(tree_in$tip.label, paste(".",gene,"..",sep=""),"")

#Make tree list + root + remove outgroup + make ultrametric
trees <- list(tree_in)
trees <- lapply(trees, function(t){root.phylo(t, outgroup = "CAR")})
trees <- lapply(trees, function(t){drop.tip(t, c(13,14))})
trees <- lapply(trees, chronos)

#Visualize
#par(mfrow = c(2,2))
#lapply(trees, plot)

#branch length distances
Dists <- lapply(trees, cophenetic.phylo)
Dists <- lapply(Dists, function(D){diag(D) <- rep(NA, length(diag(D))); return(D)})

#Get BL distances within and between populations
pops <- unique(trees[[1]]$tip.label)
within <- lapply(pops, function(P){lapply(Dists, function(D){mean(D[which(rownames(D) == P),which(colnames(D) == P)], na.rm = T)})})
between <- lapply(pops, function(P){lapply(Dists, function(D){mean(c(D[which(rownames(D) == P),which(colnames(D) != P)], D[which(rownames(D) != P),which(colnames(D) == P)]), na.rm = T)})})

within <- reverseList(within) %>% lapply(function(x){unlist(x) %>% mean()})
between <- reverseList(between) %>% lapply(function(x){unlist(x) %>% mean()})

#Determine balancing selection
BS <- lapply(as.list(1:length(within)), function(i){(within[[i]][1] - between[[i]][1])})

write.table(data.frame(geneid=gene, cophen_dist=unlist(BS)), file = args[3], sep = "\t",quote = F ,row.names = F, col.names = F )
