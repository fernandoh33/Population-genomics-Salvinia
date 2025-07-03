#!/bin/bash
#SBATCH --account=your-account
#SBATCH --time=1-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8

module load StdEnv/2023
module load angsd/0.940

REF=reference/S_molesta_cleaned_46_chr.fasta

mkdir geno.files sfs.files theta.files

#Creating genotype likelihoods files for each chromosome, for files including invaraint sites I chosen 100 random 100kb windows across the genome 
for chr in $(cat salvinia.chromosomes);
do angsd -bam bam.list -nInd 100 -minInd 10 -doSaf 1 -baq 1 -anc $REF -ref $REF -GL 2 -P 8 -minMapQ 30 -minQ 20 -out geno.files/$chr.all.samples.w.invariants.0.1 -r $chr:
done;

#Estimating sfs for each chromosome, -fold 1 is included because we do not have an ancestral reference to polarize the snps
for chr in $(cat salvinia.chromosomes);
do realSFS geno.files/$chr.all.samples.w.invariants.0.1.saf.idx -P 8 -fold 1 > sfs.files/$chr.sfs;
#Estimating diversity and neutrality statistics from the sfs, need to include the flag -fold 1 here too
realSFS saf2theta geno.files/$chr.all.samples.w.invariants.0.1.saf.idx -sfs sfs.files/$chr.sfs -outname theta.files/$chr.w.invariants.0.1 -fold 1;
#When statistics are estimated for each site, is possible to estimate statistics by windows (here 100kb) or a global statistic (one per chromosome in this case)
thetaStat do_stat theta.files/$chr.w.invariants.0.1.thetas.idx -win 100000 -step 100000 -outnames theta.files/$chr.theta.w100.s100;
done;

#Creating genotype likelihoods including only snps for population genomics
for chr in $(cat salvinia.chromosomes);
do angsd -bam bam.list -fai $REF.fai -ref $REF -nInd 100 -minInd 10 -doMajorMinor 1 -doPost 1 -doMaf 1 -doGlf 2 -out geno.files/$chr.all.samples.snps.maf0.01.minInd0.1 -gl 2 -minMapQ 30 -minQ 20 -minMaf 0.01 -SNP_pval 1e-6 -nThreads 8 -baq 1;
done;

#merge beagle files, the code is useful to filter beagle files using a list of positions, in this case 4-fold degenerate sites
#extract the header from any chromosome beagle file, for example Chr_10
zcat geno.files/Chr_10.all.samples.snps.maf0.01.minInd0.1.beagle.gz | head -1 > header.beagle
#filter the beagle files or just uncompressed them to merge
for chr in $(cat salvinia.chromosomes);
do sed 's/:/_/g' $chr.4fold.positions > tmp.$chr.4fold.positions;
zcat geno.files/$chr.all.samples.snps.maf0.01.minInd0.1.beagle.gz | awk 'NR==FNR{a[$1]; next} $1 in a' tmp.$chr.4fold.positions -> tmp.$chr.filtered.txt;
cat header.beagle tmp.$chr.filtered.txt |bgzip -c > geno.files/$chr.4fold.positions.maf.0.01.beagle.gz;
done;

#now merge all beagle files and remove temporary files
cat header.beagle tmp*.filtered.txt | bgzip -c > all.chrs.4fold.positions.maf.0.01.beagle.gz
rm tmp*

#pcangsd using 4-fold degenerate sites and three distinct maf;
for maf in 0.01 0.05 0.1;
do pcangsd -b all.chrs.4fold.positions.maf.0.01.beagle.gz -t 8 --maf $maf -o out.pcangsd.4fold.sites.maf$maf.minInd0.1;
done;
