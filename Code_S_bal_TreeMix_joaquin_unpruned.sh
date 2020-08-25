#!/bin/bash
#SBATCH -J S_bal_TreeMix
#SBATCH -t 24:00:00
#SBATCH --mem=32GB
#SBATCH -n 4

# Code for treemix

module load gsl/2.5
module load boost/1.68
module load plink/1.90
module load vcftools/0.1.16
module load bcftools/1.9
module load tabix/0.2.6

# filter sites by quality
cp ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf.gz ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k1000.vcf.gz
cp ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf.gz ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k500.vcf.gz
cp ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf.gz ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k200.vcf.gz
cp ../data/vcf_rehead/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead.vcf.gz ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k100.vcf.gz

# create treemix input files
sh vcf2treemix.sh ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k1000.vcf.gz S_bal_TreeMix_joaquin_unpruned.clust
sh vcf2treemix.sh ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k500.vcf.gz S_bal_TreeMix_joaquin_unpruned.clust
sh vcf2treemix.sh ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k200.vcf.gz S_bal_TreeMix_joaquin_unpruned.clust
sh vcf2treemix.sh ../data/vcf_split/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k100.vcf.gz S_bal_TreeMix_joaquin_unpruned.clust

# run treemix with migration edges
for i in {0..8}
do
	treemix -m $i -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k1000.treemix.frq.gz -o ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k1000.treemix.$i -root WCAN -bootstrap -k 1000 -noss  #  > treemix_$i_log
	treemix -m $i -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k500.treemix.frq.gz -o ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k500.treemix.$i -root WCAN -bootstrap -k 500 -noss  #  > treemix_$i_log
	treemix -m $i -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k200.treemix.frq.gz -o ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k200.treemix.$i -root WCAN -bootstrap -k 200 -noss  #  > treemix_$i_log
	treemix -m $i -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k100.treemix.frq.gz -o ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k100.treemix.$i -root WCAN -bootstrap -k 100 -noss  #  > treemix_$i_log
done
