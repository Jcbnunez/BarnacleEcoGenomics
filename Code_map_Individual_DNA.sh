#!/bin/bash

#########################################################

# This is an example code to map reads to a reference -- individual DNA

####################
#LOAD OSCAR MODULES#
####################

module load bwa/0.7.15 #Align reads
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
module load qualimap

#########################################################
#Script variables													<----- USER DEFINED
CPU=6
JAVAMEM=110g
QUAL=35
PRIOR=2.8e-09
#########################################################

#############################
#Input files and Identifiers# 										<----- USER DEFINED
#############################

#Input Reads 
F=.../F.trim.fastq.gz
R=.../R.trim.fastq.gz

#Name of the project: to be added to all outputs
ProjectName=<Project Name>

#References
reference=~.../reference.fa
#########################################################

########
#SCRIPT#
########
echo "Starting Process"
date

# check reads
head -n 1 $F $R

# Index reference comment out with #-# if previously indexed
#-#bwa index $reference

# map reads onto the reference with BWA-MEM
bwa mem -M -t $CPU $reference $F $R > $ProjectName.sam

# Raw mapping statistics
samtools flagstat --threads $CPU $ProjectName.sam > flagstats_raw_$ProjectName.txt

# reformat to Bam
samtools view -b -q $QUAL -f 0x0002 -F 0x0004 -F 0x0008 --threads $CPU  $ProjectName.sam > flt.$ProjectName.bam

# Sort with picard
java -Xmx$JAVAMEM -jar $PICARD SortSam I=flt.$ProjectName.bam O=flt.srt.$ProjectName.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

# Remove duplicates with picard
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates I=flt.srt.$ProjectName.bam O=flt.rmdp.srt.$ProjectName.bam M=$ProjectName.dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

# flagstats
samtools flagstat --threads $CPU flt.rmdp.srt.$ProjectName.bam > flagstats_$ProjectName.txt

# index
samtools index flt.rmdp.srt.$ProjectName.bam

# Quality Check
qualimap bamqc -bam flt.rmdp.srt.$ProjectName.bam  -outdir ./Qualimap_$ProjectName --java-mem-size=$JAVAMEM

rm $ProjectName.sam 
rm $ProjectName.bam
rm flt.$ProjectName.bam
rm flt.srt.$ProjectName.bam

#########################
# ===== Reference Call ##
#########################

bcftools mpileup --min-MQ $QUAL -Ou -f $reference -A flt.rmdp.srt.$ProjectName.bam | bcftools call -m --prior $PRIOR --skip-variants indels  > $ProjectName.allsites.vcf
bcftools filter -e'%QUAL<35 || %MAX(DP)<=10' $ProjectName.allsites.vcf > $ProjectName.filter.allsites.vcf

bgzip $ProjectName.filter.allsites.vcf
tabix $ProjectName.filter.allsites.vcf.gz

rm $ProjectName.allsites.vcf



echo "Finishing Process"
date
echo "cheers JCBN"