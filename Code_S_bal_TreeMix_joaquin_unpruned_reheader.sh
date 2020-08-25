#!/bin/bash
#SBATCH -J S_bal_TreeMix_reheader
#SBATCH -t 24:00:00
#SBATCH --mem=8GB
#SBATCH -n 4

# # # # Load required software\
module load bcftools/1.9
module load tabix/0.2.6

# # # VCF merge files
gunzip ../data/vcf_raw/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.vcf
bcftools reheader -s S_bal_TreeMix_joaquin_unpruned_sample_names.txt -o ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf ../data/vcf_raw/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.vcf
awk '{gsub(/^scaffold/,""); print}' ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf >| ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.temp.vcf
mv -f ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.temp.vcf ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf
bgzip -f ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf
tabix ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf.gz
bgzip -f ../data/vcf_raw/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.vcf
