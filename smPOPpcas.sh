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
VCF=/mnt/research/josephslab/Maya/CBP_CRCGCONP_maf_final.vcf.gz
OUTDIR=/mnt/scratch/wils1582/admixture_out
EU_SAMPLES=/mnt/home/wils1582/capsella_population_structure/vcf_neu.txt
NYC_SAMPLES=/mnt/home/wils1582/capsella_population_structure/vcf_cbp_nyc.txt
AS_SAMPLES=/mnt/home/wils1582/capsella_population_structure/vcf_easia.txt
#
# Move to directory
cd $OUTDIR

# Set up modules
module purge
ml PLINK/2.00a3.7-gfbf-2023a

# Combine the European and NYC samples
cat $EU_SAMPLES $NYC_SAMPLES > N_eurasian.txt

# Prune and PCA N. Eurasian CBP group
plink2 --vcf $VCF \
  --extract all_cbp_snps.prune.in \
  --make-bed \
  --out nyc_neu \
	--allow-extra-chr \
	--keep N_eurasian.txt \
	--set-all-var-ids @:# \
	--pca
	
# Generate frequency file for E. Asian
plink2 --vcf $VCF \
	--extract all_cbp_snps.prune.in \
	--freq 'zs' \
  --out easia_pruned_cbp \
	--allow-extra-chr \
	--keep $AS_SAMPLES \
	--set-all-var-ids @:#
	
	
# Caluclate eigenstuff for East Asian
plink2 --vcf $VCF \
	--extract all_cbp_snps.prune.in \
	--read-freq easia_pruned_cbp.afreq.zst \
  --make-bed \
  --out easia_cbp \
	--allow-extra-chr \
	--keep $AS_SAMPLES \
	--set-all-var-ids @:# \
	--pca