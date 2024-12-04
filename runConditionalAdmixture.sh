#!/bin/bash
#
#SBATCH --job-name=EUAconditionalAdmix
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-30:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A_%a.out

# Performs LD pruning on all Capsella bursa-pastoris
# Does PCA on all Capsella bursa-pastoris and C. bursa-pastoris from NYC only
# Calculates admixture proportions in Eurasia
# Calculates NYC admicture proportions conditional on Eurasian allele frequencies

### VARIABLES
VCF=/mnt/research/josephslab/Maya/capsella/vcf/filtered/CBP_CRCGCONP_maf_final.vcf.gz
OUTDIR=/mnt/scratch/wils1582/admixture_out
EU_SAMPLES=/mnt/home/wils1582/capsella_population_structure/vcf_cbp_eurasia.txt
NYC_SAMPLES=/mnt/home/wils1582/capsella_population_structure/vcf_cbp_nyc.txt
ALL_CBP=/mnt/home/wils1582/capsella_population_structure/all_cbp.txt

# move to directory
#mkdir -p $OUTDIR
cd $OUTDIR

#load modules
module purge
ml PLINK/2.00a3.7-gfbf-2023a ADMIXTURE/1.3.0
#
##### LD Prune and output BED
## Identify SNPs in LD
## Parameters consider SNPs in 250kb windows, cutoff is and r^2 correlation of 0.125 between a pair of SNPs. Step size is 1 because PLINK2.0 says so
## since I specified the window size in kilobases

# Now doing 500 snps with 50 step size, and r^2 of 0.5 from Shirsekar et. al. 2021
plink2 --vcf $VCF \
	--indep-pairwise 500 50 0.5 \
	--keep $ALL_CBP \
	--allow-extra-chr \
	--out all_cbp_snps \
	--set-all-var-ids @:# 

# Prune LD SNPs
# Variant IDs set to be 'scaffold:position'
# ADMIXTURE requires the input be in PLINK .bed format, so here I create the bed file as well.
# LD pruning anf PCA of all CBP together
plink2 --vcf $VCF \
	--extract all_cbp_snps.prune.in \
  --make-bed \
  --out all_pruned_cbp \
        --allow-extra-chr \
        --keep $ALL_CBP \
        --set-all-var-ids @:# \
        --pca
# Prune Eurasian CBP
plink2 --vcf $VCF \
	--extract all_cbp_snps.prune.in \
  --make-bed \
  --out eurasian_pruned_cbp \
	--allow-extra-chr \
	--keep $EU_SAMPLES \
	--set-all-var-ids @:# \
	--pca
# Prune New York/New Jersey CBP
plink2 --vcf $VCF \
  --extract all_cbp_snps.prune.in \
  --make-bed \
  --out nyc_pruned_cbp \
	--allow-extra-chr \
	--keep $NYC_SAMPLES \
	--set-all-var-ids @:# \
	--pca
# Check if the .bim files are the same
diff -s eurasian_pruned_cbp.bim nyc_pruned_cbp.bim
# ADMIXTURE only takes chromosomes as numbers so in all the PLINK output files, jlSCF_ needs to be removed
sed -i 's/^jlSCF_//g' eurasian_pruned_cbp.bed
sed -i 's/^jlSCF_//g' eurasian_pruned_cbp.bim
sed -i 's/^jlSCF_//g' eurasian_pruned_cbp.fam
sed -i 's/^jlSCF_//g' nyc_pruned_cbp.bed
sed -i 's/^jlSCF_//g' nyc_pruned_cbp.bim
sed -i 's/^jlSCF_//g' nyc_pruned_cbp.fam

# Also change the names of the final two contigs to be numeric only; makes them called 21 and 22
# sed -i 's/^contig_/2/g' eurasian_pruned_cbp.bed
# sed -i 's/^contig_/2/g'	eurasian_pruned_cbp.bim
# sed -i 's/^contig_/2/g'	eurasian_pruned_cbp.fam
# sed -i 's/^contig_/2/g' nyc_pruned_cbp.bed
# sed -i 's/^contig_/2/g' nyc_pruned_cbp.bim
# sed -i 's/^contig_/2/g' nyc_pruned_cbp.fam

# Run unsupervised ADMIXTURE with K=2
admixture --cv eurasian_pruned_cbp.bed 2 | tee eua_log2.out

# Use learned allele frequencies as (fixed) input to next step
cp eurasian_pruned_cbp.2.P nyc_pruned_cbp.2.P.in

# Run projection ADMIXTURE with K=2
admixture -P nyc_pruned_cbp.bed 2

### K = 3
# Run unsupervised ADMIXTURE with K=3
admixture --cv eurasian_pruned_cbp.bed 3 | tee eua_log3.out

# Use learned allele frequencies as (fixed) input to next step
cp eurasian_pruned_cbp.3.P nyc_pruned_cbp.3.P.in

# Run projection ADMIXTURE with K=3
admixture -P nyc_pruned_cbp.bed 3

### K = 4
# Run unsupervised ADMIXTURE with K=4
admixture --cv eurasian_pruned_cbp.bed 4 | tee eua_log4.out

# Use learned allele frequencies as (fixed) input to next step
cp eurasian_pruned_cbp.4.P nyc_pruned_cbp.4.P.in

# Run projection ADMIXTURE with K=4
admixture -P nyc_pruned_cbp.bed 4

### K = 5
# Run unsupervised ADMIXTURE with K=5
admixture --cv eurasian_pruned_cbp.bed 5 | tee eua_log5.out

# Use learned allele frequencies as (fixed) input to next step
cp eurasian_pruned_cbp.5.P nyc_pruned_cbp.5.P.in

# Run projection ADMIXTURE with K=5
admixture -P nyc_pruned_cbp.bed 5

### K = 6
# Run unsupervised ADMIXTURE with K=6
admixture --cv eurasian_pruned_cbp.bed 6 | tee eua_log6.out

# Use learned allele frequencies as (fixed) input to next step
cp eurasian_pruned_cbp.6.P nyc_pruned_cbp.6.P.in

# Run projection ADMIXTURE with K=6
admixture -P nyc_pruned_cbp.bed 6


# write cv errors to file
grep -h CV eua_log*.out > eua_cv_error.log
