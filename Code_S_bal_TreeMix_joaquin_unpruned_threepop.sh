#!/bin/bash
#SBATCH -J S_bal_TreeMix
#SBATCH -t 24:00:00
#SBATCH --mem=32GB
#SBATCH -n 4

module load gsl/2.5
module load boost/1.68
module load plink/1.90
module load vcftools/0.1.16
module load bcftools/1.9
module load tabix/0.2.6

threepop -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k1000.treemix.frq.gz -k 1000 >| ../results/threepop/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k1000.treemix.threepop.txt
threepop -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k500.treemix.frq.gz -k 500 >| ../results/threepop/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k500.treemix.threepop.txt
threepop -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k200.treemix.frq.gz -k 200 >| ../results/threepop/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k200.treemix.threepop.txt
threepop -i ../results/treemix/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k100.treemix.frq.gz -k 100 >| ../results/threepop/All_pops_Sbal_VCF.filter.99misd.maf5.SNPs.rehead_k100.treemix.threepop.txt
