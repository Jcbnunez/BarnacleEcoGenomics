#!/bin/bash

#SBATCH -J Zonatated_Theta_ME
#SBATCH -n 1
#SBATCH --mem=90GB
#SBATCH -t 90:00:00

################
# Load Modules #
################

module load samtools/1.9
module load R/3.4.3

###################
# Comp. Variables #
###################

CPU=1
JAVAMEM=80g

#####################
# Pipeline programs # ######## ####### ######## do not change (unless for debug)
#####################

Identify_indels=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/basic-pipeline/identify-genomic-indel-regions.pl
Filter_pileup_by_GTF=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/basic-pipeline/filter-pileup-by-gtf.pl
Subsample_pileup=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/basic-pipeline/subsample-pileup.pl
Variance_sliding=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/Variance-sliding.pl

###################
# External files  #
###################

Panel=~/data/S_balanoides_genomics_resources/Analyses/Zonation_Paper/0.Merged_bams_panel/Hodgsons_Panel.bam

#Name of out files
POP=ME
pool_sizes=288

#####################
# Script. Variables #
#####################

ProjectName=Zonation_est

reference=~/data/S_balanoides_genomics_resources/Genomes/25XTENXME_10XPACBIO/DBG2OLC_Sparce_SGAcorrectedANDFiltered_GOOD/Sbal3_redundans/redundans/scaffolds.reduced.fa

#parameters
WS=50
SS=10
MAXCOV=600
target_cov=100

ZONATED_regions=~/data/S_balanoides_genomics_resources/Analyses/Zonation_Paper/7.Zonated_Theta_est/ZONATION_regions.bed

#####################
# Pipeline programs # ######## ####### ######## do not change (unless for debug)
#####################

Identify_indels=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/basic-pipeline/identify-genomic-indel-regions.pl
Filter_pileup_by_GTF=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/basic-pipeline/filter-pileup-by-gtf.pl
Subsample_pileup=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/basic-pipeline/subsample-pileup.pl
Variance_sliding=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/Variance-sliding.pl
Var_at_Pos=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation1/Variance-at-position.pl
 
#################################

# Part 0: Create log file
touch $ProjectName.$POP.log
echo $(date) "Starting process ${ProjectName}.${POP}" >> $ProjectName.$POP.log

# Start the step counter
let i=1

#################################

# Part 1: Make pileup file
 echo $(date)  "Starting step ${i}: making a mpileup of the desired populations" >> $ProjectName.$POP.log
samtools mpileup -q 35 -Q 20 -f $reference -l $ZONATED_regions $Panel > $ProjectName.$POP.pile

#Exit status check
	if [ -s $ProjectName.$POP.pile ]
		then
		echo $(date) "Step ${i} seems ok, moving on"  >> $ProjectName.$POP.log
		((++i))
		
		else	
		echo $(date) "the file is empty! Exiting pipe now!" >> $ProjectName.$POP.log
		echo $(date) "check the script, debug, and rerun" >> $ProjectName.$POP.log
		echo $(date) "Pipeline died at step ${i}" >> $ProjectName.$POP.log
	    exit 1
	fi


#################################

# Part 2: Identify indels
echo $(date)  "Starting step ${i}: Identify indels" >> $ProjectName.$POP.log

perl $Identify_indels --input $ProjectName.$POP.pile --output $ProjectName.$POP.indel.regions.gtf --indel-window 5

#Exit status check
	if [ -s $ProjectName.$POP.indel.regions.gtf ]
		then
		echo $(date) "Step ${i} seems ok, moving on"  >> $ProjectName.$POP.log
		((++i))
		
		else	
		echo $(date) "the file is empty! Exiting pipe now!" >> $ProjectName.$POP.log
		echo $(date) "check the script, debug, and rerun" >> $ProjectName.$POP.log
		echo $(date) "Pipeline died at step ${i}" >> $ProjectName.$POP.log
	    exit 1
	fi

#################################

 # Part 3: Removing indels
echo $(date)  "Starting step ${i}: Remove indels" >> $ProjectName.$POP.log

perl $Filter_pileup_by_GTF --gtf  $ProjectName.$POP.indel.regions.gtf --input $ProjectName.$POP.pile --output $ProjectName.$POP.noindels.pile

#Exit status check
	if [ -s $ProjectName.$POP.noindels.pile ]
		then
		echo $(date) "Step ${i} seems ok, moving on"  >> $ProjectName.$POP.log
		((++i))
		
		else	
		echo $(date) "the file is empty! Exiting pipe now!" >> $ProjectName.$POP.log
		echo $(date) "check the script, debug, and rerun" >> $ProjectName.$POP.log
		echo $(date) "Pipeline died at step ${i}" >> $ProjectName.$POP.log
	    exit 1
	fi

#################################

# Part 4: Subsample coverage
echo $(date)  "Starting step ${i}: subsample coverage" >> $ProjectName.$POP.log

perl  $Subsample_pileup --min-qual 20 --method withoutreplace --max-coverage $MAXCOV --fastq-type sanger --target-coverage $target_cov --input $ProjectName.$POP.noindels.pile --output $ProjectName.$POP.Cov_$target_cov.noindels.pile

 #Exit status check
	if [ -s $ProjectName.$POP.Cov_$target_cov.noindels.pile ]
		then
		echo $(date) "Step ${i} seems ok, moving on"  >> $ProjectName.$POP.log
		((++i))
		
		else	
		echo $(date) "the file is empty! Exiting pipe now!" >> $ProjectName.$POP.log
		echo $(date) "check the script, debug, and rerun" >> $ProjectName.$POP.log
		echo $(date) "Pipeline died at step ${i}" >> $ProjectName.$POP.log
	    exit 1
	fi

#################################

# Part 5: Estimate pi
echo $(date)  "Starting step ${i}: estimate pi" >> $ProjectName.$POP.log


perl $Variance_sliding --input $ProjectName.$POP.Cov_$target_cov.noindels.pile  --output $ProjectName.$POP.pi.$WS.$SS.txt --measure pi --window-size $WS --step-size $SS --min-count 2 --min-coverage 10 --max-coverage 100 --pool-size $pool_sizes --min-qual 20 --fastq-type sanger

  #Exit status check
	if [ -s $ProjectName.$POP.pi.$WS.$SS.txt ]
		then
		echo $(date) "Step ${i} seems ok, moving on"  >> $ProjectName.$POP.log
		((++i))
		
		else	
		echo $(date) "the file is empty! Exiting pipe now!" >> $ProjectName.$POP.log
		echo $(date) "check the script, debug, and rerun" >> $ProjectName.$POP.log
		echo $(date) "Pipeline died at step ${i}" >> $ProjectName.$POP.log
	    exit 1
	fi

#################################

# Part 6: Estimate D

echo $(date)  "Starting step ${i}: estimate D" >> $ProjectName.$POP.log

perl $Variance_sliding --input $ProjectName.$POP.Cov_$target_cov.noindels.pile  --output $ProjectName.$POP.D.$WS.$SS.txt --measure D --window-size $WS --step-size $SS --min-count 2 --min-coverage 10 --max-coverage 100 --pool-size $pool_sizes --min-qual 20 --fastq-type sanger


  #Exit status check
	if [ -s $ProjectName.$POP.D.$WS.$SS.txt ]
		then
		echo $(date) "Step ${i} seems ok, I am done!"  >> $ProjectName.$POP.log
		
		else	
		echo $(date) "the file is empty! Exiting pipe now!" >> $ProjectName.$POP.log
		echo $(date) "check the script, debug, and rerun" >> $ProjectName.$POP.log
		echo $(date) "Pipeline died at step ${i}" >> $ProjectName.$POP.log
	    exit 1
	fi

#################################

echo "Finished Loop"
date