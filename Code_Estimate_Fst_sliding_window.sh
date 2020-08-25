#!/bin/bash

#SBATCH -J Make_Fst_allelFreq_Genes
#SBATCH -n 1
#SBATCH --mem=60GB
#SBATCH -t 120:00:00

date

#Load Modules
module load samtools/1.9

# outputfile identifier <<<+
Project=Sbal_habWide_MKT

##################
# IN- FILES ######
##################

SYNC=~/data/S_balanoides_genomics_resources/Analyses/Genome_Paper/14.fst/Sbal_habWide_MKT.noindels.nosing.fstpops.sync

##########################
# Uder defined variables # ##### change according to input
##########################

GUIDE=~/data/S_balanoides_genomics_resources/Gene_prediction_models/Fly_RNA_prior/chr_genes_bed.bed
reference=~/data/S_balanoides_genomics_resources/Genomes/25XTENXME_10XPACBIO/DBG2OLC_Sparce_SGAcorrectedANDFiltered_GOOD/Sbal3_redundans/redundans/scaffolds.reduced.fa

JAVAMEM=50g
WS=5000
SS=2500

#####################
# Pipeline programs # ######## ####### ######## do not change (unless for debug)
#####################

fst_sliding=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation2/fst-sliding.pl
snp_frequency_diff=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation2/snp-frequency-diff.pl

###################
# Run the program # ######## ####### ######## do not change (unless for debug)
###################

##### ######
# RUN PIPE #
#### ##### #

#order pops
# ME RI ICE NOR UKW WCAN
#37:38:20:20:28:20 (hap size)
#74:76:40:40:56:40 (dip size)

#perl $fst_sliding --input $SYNC --output $Project.$WS.$SS.fst.txt --min-count 6 --min-coverage 10 --max-coverage 100 --window-size $WS --step-size $SS --pool-size 74:76:40:40:56:40

perl $snp_frequency_diff --input $SYNC --output $Project.Exons.AlleleCount.txt --min-count 5 --min-coverage 10 --max-coverage 200


date

##############
## END PIPE ##
##############