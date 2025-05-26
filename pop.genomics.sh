#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

REF=reference/S_molesta_cleaned_46_chr.fasta

mkdir geno.files sfs.files theta.files

#Creating genotype likelihoods files for each chromosome, for files including invaraint sites I generate one per chromosome only because file are too large
for chr in $(cat salvinia.chromosomes);
do angsd -bam bam.list -nInd 100 -minInd 10 -doSaf 1 -baq 1 -anc $REF -ref $REF -GL 2 -P 8 -minMapQ 30 -minQ 20 -out geno.files/$chr.all.samples.w.invariants.0.1 -r $chr:
done;

#Estimating sfs for each chromosome, make sure to include -fold 1 if you don't have an ancestral reference to polarize the snps
for chr in $(cat salvinia.chromosomes);
do realSFS geno.files/$chr.all.samples.w.invariants.0.1.saf.idx -P 8 -fold 1 > sfs.files/$chr.sfs;
#Estimating diversity and neutrality statistics from the sfs, include -fold 1 here too
realSFS saf2theta geno.files/$chr.all.samples.w.invariants.0.1.saf.idx -sfs sfs.files/$chr.sfs -outname theta.files/$chr.w.invariants.0.1 -fold 1;
#When statistics are estiamated for each site, is possible to estimate statistics by windows (here 100kb) or a global statistic (one per chromosome in this case)
thetaStat do_stat theta.files/$chr.w.invariants.0.1.thetas.idx -win 100000 -step 100000 -outnames theta.files/$chr.theta.w100.s100;
done;

#Creating genotype likelihoods including only snps for population genomics, these files can be very large too
angsd -bam bam.list -fai $REF.fai -ref $REF -nInd 100 -minInd 10 -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -out all.samples.snps.maf0.01.minInd0.1 -gl 2 -minMapQ 30 -minQ 20 -minMaf 0.01 -SNP_pval 1e-6 -nThreads 8 -baq 1;

#Ideally snps should be filtered by low LD for population genomic analyses, though in this case all samples seems to be a single clone with very low diversity, so we are using all snps
#Run pcangsd with the snps dataset, testing different maf filtering is recommended;
for maf in 0.01 0.05 0.1;
do pcangsd -b all.samples.snps.maf0.01.minInd0.1.beagle.gz -t 8 --maf $maf -o out.pcangsd.maf$maf.minInd0.1;
done

#Run NGSadmix
for i in `seq 1 5`;do NGSadmix -likes all.samples.snps.maf0.01.minInd0.1.beagle.gz -K $i -outfiles out.admixure.maf0.01.minInd0.1.K$i -P 12;done
