#!/bin/bash

module load samtools
module load vcftools
module load bcftools
module load tabix/0.2.6 
module load iq-tree/1.6.7 
module load R/3.6.0  ## <<< All packages function in this version >>>> ####
module load gcc
module load imagemagick/7.0.7


GUIDE=.../regions_of_interest.bed
reference=.../reference_genome.fasta
QUAL=<QUAL>
BO=<Bootstrap N for trees>

touch final.cophen.out.txt
mkdir trees
mkdir fasta_files
mkdir btrees

echo "start"

awk -F "\t"  '{print $1,"\t",$2}' $GUIDE | while read POS GENE
 do 
 echo ${POS} "is" ${GENE}
 
 #############################
 #Canada
 #############################
 
 echo "Canada"
 
 BAM=""
 VCF=""
 POP=""
 
 BAM=~/...bam
 VCF=~/...vcf.gz
 POP="PopName" # ... in CAR CAN NOR UKW ICE ME RI
 
 ####### PROGRAMS #####
 HAPCUT=.../HAPCUT/HapCUT2/build
 
 ####### CODE #####
 samtools view -b -h $BAM ${POS} > ${POP}.${GENE}.${POS}.bam
 
 bcftools mpileup --min-MQ $QUAL -Ou -f $reference -A ${POP}.${GENE}.${POS}.bam | bcftools call -cv --skip-variants indels  > ${POP}.${GENE}.${POS}.vcf
 
 $HAPCUT/extractHAIRS --bam ${POP}.${GENE}.${POS}.bam --VCF ${POP}.${GENE}.${POS}.vcf --out ${POP}.${GENE}.${POS}.fragment_file
 
 $HAPCUT/HAPCUT2 --fragments ${POP}.${GENE}.${POS}.fragment_file --VCF ${POP}.${GENE}.${POS}.vcf --output ${POP}.${GENE}.${POS}.haplotype_output_file
 
 sed 's:1/1:1|1:g' ${POP}.${GENE}.${POS}.haplotype_output_file.phased.VCF > ${POP}.${GENE}.${POS}.haplotype_output_file.1.phased.VCF
 
 
 vcftools --vcf ${POP}.${GENE}.${POS}.haplotype_output_file.1.phased.VCF --min-alleles 2 --max-alleles 2 --minQ $QUAL --minDP 6  --phased  --recode --recode-INFO-all --out ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly
 
 
 bgzip ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly.recode.vcf
 tabix ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly.recode.vcf.gz
 
 
 bcftools consensus -H 1 -f <(samtools faidx $reference ${POS}) ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly.recode.vcf.gz > ${POP}.${GENE}.${POS}.haplotype0.fasta
  
 bcftools consensus -H 2 -f <(samtools faidx $reference ${POS}) ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly.recode.vcf.gz > ${POP}.${GENE}.${POS}.haplotype1.fasta
 
 sed -e 's/^>'${POS}'/>'${POP}.${GENE}'.0/g' ${POP}.${GENE}.${POS}.haplotype0.fasta > ${POP}.${GENE}.${POS}.haplotype0.named.fasta
 
 sed -e 's/^>'${POS}'/>'${POP}.${GENE}'.1/g' ${POP}.${GENE}.${POS}.haplotype1.fasta > ${POP}.${GENE}.${POS}.haplotype1.named.fasta
 
 rm ${POP}.${GENE}.${POS}.bam
 rm ${POP}.${GENE}.${POS}.vcf
 rm ${POP}.${GENE}.${POS}.fragment_file
 rm ${POP}.${GENE}.${POS}.haplotype_output_file
 rm ${POP}.${GENE}.${POS}.haplotype0.fasta
 rm ${POP}.${GENE}.${POS}.haplotype1.fasta
 rm ${POP}.${GENE}.${POS}.haplotype_output_file.1.phased.VCF
 rm ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly.recode.vcf.gz.tbi
 rm ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly.recode.vcf.gz
 rm ${POP}.${GENE}.${POS}.haplotype_output_file.phasedOnly.log
 rm ${POP}.${GENE}.${POS}.haplotype_output_file.phased.VCF
 
 
#############################
 #Merge all
 #############################
 
 cat *.fasta > ${GENE}.fasta
 
 for POP in CAR CAN NOR UKW ICE ME RI
 do
 rm ${POP}.${GENE}.${POS}.haplotype0.named.fasta
 rm ${POP}.${GENE}.${POS}.haplotype1.named.fasta
 done
 
 #############################
 #Build tree
 #############################

cp ./fasta_files/${GENE}.fasta ./

iqtree -s ${GENE}.fasta -bo $BO -o CAR.${GENE}.0 -nt AUTO

rm ${GENE}.fasta.treefile
rm ${GENE}.fasta.splits.nex
rm ${GENE}.fasta.log
rm ${GENE}.fasta.iqtree
rm ${GENE}.fasta.ckp.gz
rm ${GENE}.fasta.model.gz
rm ${GENE}.fasta.mldist
rm ${GENE}.fasta.bionj
rm ${GENE}.fasta.contree
rm ${GENE}.fasta

mv ${GENE}.fasta.boottrees ./btrees
