#!/bin/bash

#SBATCH -J Make_Fst_allelFreq_Genes
#SBATCH -n 1
#SBATCH --mem=60GB
#SBATCH -t 120:00:00

date

#Load Modules
module load samtools/1.9

# outputfile identifier <<<+
Project=<Name of Project>

##################
# IN- FILES ######
##################

SYNC=~/popoolation/sync

##########################
# Uder defined variables # ##### change according to input
##########################

GUIDE=~/.../chr_genes_bed.bed #Chromosomes to investigate
reference=/.../refererence.fasta

JAVAMEM=50g
WS=<Window Size>
SS=<Step Size>
POP_SIZES=<Pop:Sizes:of:all:pops:...>

#####################
# Pipeline programs # ######## ####### ######## do not change (unless for debug)
#####################

fst_sliding=.../popoolation2/fst-sliding.pl
snp_frequency_diff=.../popoolation2/snp-frequency-diff.pl

###################
# Run the program # ######## ####### ######## do not change 
###################

##### ######
# RUN PIPE #
#### ##### #

perl $fst_sliding --input $SYNC --output $Project.$WS.$SS.fst.txt --min-count 6 --min-coverage 10 --max-coverage 100 --window-size $WS --step-size $SS --pool-size $POP_SIZES

perl $snp_frequency_diff --input $SYNC --output $Project.Exons.AlleleCount.txt --min-count 5 --min-coverage 10 --max-coverage 200

date

##############
## END PIPE ##
##############