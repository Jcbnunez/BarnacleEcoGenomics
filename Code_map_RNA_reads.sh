#!/bin/bash

#SBATCH -J HISAT2_HighD_RNAexp
#SBATCH -n 4
#SBATCH --mem=60GB
#SBATCH -t 30:00:00

#########################################################

# This is an example code to map reads to a reference -- individual RNA

##############
#SYSTEM CHECK#
##############
echo "-----------------------------"
echo ">>>>The current folder is<<<<"
pwd
echo ">>>>Operating System Information<<<<"
uname -a
echo ">>>>I/O Information<<<<"
ulimit -a
echo ">>>>Computing Node is<<<<"
hostname
echo ">>>>Domain Name is<<<<"
domainname
echo ">>>>System DNS is<<<<"
dnsdomainname
echo "-----------------------------"

#########################################################

####################
#LOAD OSCAR MODULES#
####################

module load hisat2/2.1.0 #Align reads
module load samtools/1.9 # Manipulate SAM/BAM files
module load bcftools/1.9
module load picard-tools/2.9.2  #Sort and remove duplicates
module load bedtools/2.26.0 
module load bbmap/38.23
module load fgbio/0.6.1 
module load tabix/0.2.6 
module load vcftools/0.1.16
module load seqtk/1.3
module load gatk/4.0.9.0 

#Additional features & Programs
qualimap=~/data/Jcbn/Software/qualimap_v2.1/qualimap
HAPCUT=/users/jnunez/software/HapCUT2/build

#########################################################
#Script variables
CPU=4
JAVAMEM=110g
QUAL=20
PRIOR=2.8e-09
#########################################################

#########################
#DECLARE INPUT VARIABLES#
#########################

#Input Reads 
F=.../F.trim.fastq.gz
R=.../R.trim.fastq.gz

#Name of the project: to be added to all outputs
ProjectName=<Project Name>

#References
reference=~.../reference.fa

#########################################################

ref_pos=~/data/S_balanoides_genomics_resources/bam_files/single_individuals_RNA/

########
#SCRIPT#
########
echo "Starting Process ..."
date

# call SNPs with sam tools
# the --dta flag is mandatory for differential expression analysis
hisat2 -x $ref_pos/Sbal3_ref -1 $F -2 $R -p $CPU --dta -S $ProjectName.sam

# Raw mapping statistics
samtools flagstat --threads $CPU $ProjectName.sam > flagstats_raw_$ProjectName.txt

# reformat to Bam
samtools view -q $QUAL -F 0x0004 --threads $CPU -b $ProjectName.sam > flt.$ProjectName.bam

# Remove non-mapping reads, fix paired SAM flags added by BWA
samtools flagstat --threads $CPU flt.$ProjectName.bam > flagstats_filter_$ProjectName.txt

# Sort with picard(DNA/RNA)
java -Xmx$JAVAMEM -jar $PICARD SortSam I=flt.$ProjectName.bam O=flt.srt.$ProjectName.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

#flagstats
samtools flagstat --threads $CPU flt.srt.$ProjectName.bam > flagstats_final_$ProjectName.txt

#index
samtools index flt.srt.$ProjectName.bam

#Quality Check
module load qualimap/2.2.1 

qualimap bamqc -bam flt.srt.$ProjectName.bam  -outdir ./Qualimap_$ProjectName --java-mem-size=$JAVAMEM

rm $ProjectName.sam 
rm $ProjectName.bam
rm flt.$ProjectName.bam

#########################################################

###########################
#Begin SNP calling to VCF #
###########################
reference=~/data/S_balanoides_genomics_resources/Genomes/25XTENXME_10XPACBIO/DBG2OLC_Sparce_SGAcorrectedANDFiltered_GOOD/Sbal3_redundans/redundans/scaffolds.reduced.fa

bcftools mpileup --min-MQ $QUAL -Ou -f $reference -A flt.srt.$ProjectName.bam | bcftools call -mv --prior $PRIOR --skip-variants indels  > $ProjectName.vcf
bgzip $ProjectName.vcf
tabix $ProjectName.vcf.gz

echo "Finishing Process"
date
echo "cheers JCBN"
