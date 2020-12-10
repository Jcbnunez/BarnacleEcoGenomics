#!/bin/bash

#SBATCH -J Make_SYNC_file
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

CPU=3
JAVAMEM=80g

popoolation2=~/data/S_balanoides_genomics_resources/Misc_resources/popoolation2

###################
# External files  #
###################

reference=~/data/S_balanoides_genomics_resources/Genomes/25XTENXME_10XPACBIO/DBG2OLC_Sparce_SGAcorrectedANDFiltered_GOOD/Sbal3_redundans/redundans/scaffolds.reduced.fa

#Maine Pools

HILC=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Hodson_Island_Low_Cold/flt.rmdp.srt.ME_HI_LC.sbal3.1.bam
HIHC=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Hodson_Island_High_Cold/flt.rmdp.srt.ME_HI_HC.sbal3.1.bam
HILH=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Hodson_Island_Low_Hot/flt.rmdp.srt.ME_HI_LH.sbal3.1.bam
HIHH=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Hodson_Island_High_Hot/flt.rmdp.srt.ME_HI_HH.sbal3.1.bam

FILC=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Farmers_Island_Low_Cold/flt.rmdp.srt.ME_FI_LC.sbal3.1.bam
FIHC=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Farmers_Island_High_Cold/flt.rmdp.srt.ME_FI_HC.sbal3.1.bam
FILH=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Farmers_Island_Low_Hot/flt.rmdp.srt.ME_FI_LH.sbal3.1.bam
FIHH=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Maine_Damariscotta_High_Low_tide/Farmers_Island_High_Hot/flt.rmdp.srt.ME_FI_HH.sbal3.1.bam

#Rhode Island Pool

RICLN=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Rhode_Island_Jamestown_High_Low_tide/Site_C_Low_North/flt.rmdp.srt.RI_C_LC.sbal3.1.bam
RICHS=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Rhode_Island_Jamestown_High_Low_tide/Site_C_High_South/flt.rmdp.srt.RI_C_HH.sbal3.1.bam

RIDLN=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Rhode_Island_Jamestown_High_Low_tide/Site_D_Low_North/flt.rmdp.srt.RI_D_LC.sbal3.1.bam
RIDHS=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Rhode_Island_Jamestown_High_Low_tide/Site_D_High_South/flt.rmdp.srt.RI_D_HH.sbal3.1.bam


# Norway Pools

GBLN=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Bergen_Norway/GBLN/flt.rmdp.srt.GBLN_pool.sbal3.1.bam
GBHS=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Bergen_Norway/GBHS/flt.rmdp.srt.GBHS_pool.sbal3.1.bam

NRLN=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Bergen_Norway/NRLN/flt.rmdp.srt.NRLN_pool.sbal3.1.bam
NRHS=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Bergen_Norway/NRHS/flt.rmdp.srt.NRHS_pool.sbal3.1.bam


# Sweden Pools

YVLN=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Tjarno_Sweden/YVLN/flt.rmdp.srt.YVLN_pool.sbal3.1.bam
YVHS=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Tjarno_Sweden/YVHS/flt.rmdp.srt.YVHS_pool.sbal3.1.bam

# North Atlantic pools

ICE_pool=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Iceland_Reykjavik/flt.rmdp.srt.ICE_pool.sbal3.1.bam
UKW_pool=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Wales_UK/flt.rmdp.srt.UKW_pool.sbal3.1.bam
NOR_pool=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Norddal_Norway/flt.rmdp.srt.NOR_pool.sbal3.1.bam

# Outgroup

CAN_pool=~/data/S_balanoides_genomics_resources/bam_files/pool_seq_DNA/Canada_British_Columbia/flt.rmdp.srt.WCAN_pool.sbal3.1.bam

#####################
# Script. Variables #
#####################

# Pool sizes are as such 2*{p}; where p is ploidy = pop1:pop2:pop3:...:popN
POPSIZES=72:72:72:72:72:72:72:72:72:72:72:72:40:68:42:40

#Name of out files
Project=SCGP_Sb3.1

#parameters
#targetCOV=30
WINSIZE=1
WINSTEP=1

#Number of populations minus the out group. i.e. if 10 pop used and one is the out-group this number is 9. If no out-group was used this parameter is irreleventant. Also, if there is no out-group remember to comment out the out-group N removal portion down below.

NPOP_1=15


################################
# Tests of Popoolation modules #
################################

perl  $popoolation2/indel_filtering/identify-indel-regions.pl --test
perl  $popoolation2/indel_filtering/filter-sync-by-gtf.pl --help
java -ea -Xmx$JAVAMEM -jar $popoolation2/mpileup2sync.jar --help
perl $popoolation2/snp-frequency-diff.pl --test
perl $popoolation2/cmh-test.pl  --test

###################
# Script	      # ================> Synchronized Pipelines
###################

# Create a pileup file
samtools mpileup -f $reference -B $HIHH $FIHH $HIHC $FIHC $HILH $FILH $HILC $FILC $RICLN $RICHS $RIDHS $RIDLN $GBHS $GBLN $NRHS $NRLN $YVHS $YVLN $ICE_pool $UKW_pool $NOR_pool $CAN_pool > $Project.mpileup

#filter indels part1: create a map file
perl  $popoolation2/indel_filtering/identify-indel-regions.pl --indel-window 5 --min-count 2 --input $Project.mpileup --output $Project.gtf

#filter indels part2: filter index via gft map
perl  $popoolation2/indel_filtering/filter-sync-by-gtf.pl --input $Project.mpileup --gtf $Project.gtf --output $Project.noindels.mpileup

# Make synchronized pileup
java -ea -Xmx$JAVAMEM -jar $popoolation2/mpileup2sync.jar --input $Project.noindels.mpileup --output $Project.noindels.sync --fastq-type sanger --min-qual 35 --threads $CPU

# Eliminate Singleton calls in the sync file
# Its impossible to differentiate true singeltons/doubletons from sequencing errors. Thus, I will remove them from the Sync File.
cat $Project.noindels.sync | sed -E 's|:1:|:0:|g' | sed -E 's|\t1:|\t0:|g' | sed -E 's|:2:|:0:|g' | sed -E 's|\t2:|\t0:|g' > $Project.noindels.NoSigt.sync

#Subsample to uniform covertage
# perl ~/data/Jcbn/Software/popoolation2/subsample-synchronized.pl --input $Project.noindels.sync --output $Project.noindels.$targetCOV.sync --target-coverage $targetCOV  --max-coverage 100 --method withoutreplace
