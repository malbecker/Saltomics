#!/bin/bash
#SBATCH -N 1
#SBATCH -p macmanes,shared
#SBATCH --cpus-per-task=24
#SBATCH --job-name=cinerea_ANGSD
#SBATCH --output=cinerea_ANGSD_log.txt

#source activate orp_v2

WD=$(pwd)
TRANSCRIPTOME=$"Hyla_cinerea_transcriptome.fasta"
#PICARD=$(which picard)
ANGSD=$"/mnt/lustre/macmaneslab/ams1236/software/angsd/angsd"
cd reads/
SAMPLES=$(ls *R1_001.fastq | sed "s/R1_001.fastq//g" | grep -v subsamp)
cd ..

echo Working with the transcriptome: $TRANSCRIPTOME
echo Working with the following samples: $SAMPLES

#mkdir rcorr/

#for i in $SAMPLES
#do
#trimmomatic PE -threads 24 -baseout ${WD}/rcorr/${i}.TRIM.fastq ${WD}/reads/${i}R1_001.fastq ${WD}/reads/${i}R2_001.fastq LEADING:3 TRAILING:3 ILLUMINACLIP:/mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/barcodes/barcodes.fa:2:30:10 MINLEN:25
#done

#for i in $SAMPLES
#do
#perl /mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/software/anaconda/install/envs/orp_v2/bin/run_rcorrector.pl -t 24 -k 31 -p ${WD}/rcorr/${i}.TRIM_1P.fastq ${WD}/rcorr/${i}.TRIM_2P.fastq -od ${WD}/rcorr
#done

#mkdir bamfiles/
#cp $TRANSCRIPTOME bamfiles/

# ANGSD issue caused by age of transcriptome index, sleep to solve
#sleep 150s

#source deactivate orp_v2
module load linuxbrew/colsa

# create an BWA index from the reference transcriptome
cd bamfiles/
#bwa index $TRANSCRIPTOME

# Map individual samples to reference transcriptome using BWA
#for i in $SAMPLES; do
#    bwa mem $TRANSCRIPTOME  -t 24 \
#    ${WD}/rcorr/${i}.TRIM_1P.cor.fq \
#    ${WD}/rcorr/${i}.TRIM_1P.cor.fq > ${WD}/bamfiles/${i}.sam
#done


# run picard to produce sorted sam files
#parallel -j 4 picard SortSam INPUT={}.sam OUTPUT={}_sorted_reads.bam SORT_ORDER=coordinate ::: $SAMPLES

# Add group information to files
#sorted=$(ls *_sorted_reads.bam | sed "s/_sorted_reads.bam//g")
#parallel -j 4 picard AddOrReplaceReadGroups  I={}_sorted_reads.bam O={}_marked_groups.bam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={} ::: $SAMPLES

# Mark duplicates here
#marked=$(ls *_marked_groups.bam | sed "s/_marked_groups.bam//g")
#parallel -j 4 picard MarkDuplicates INPUT={}_marked_groups.bam OUTPUT={}_dedup_reads.bam METRICS_FILE={}_metrics.txt REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 ::: $marked

# sample list file for angsd
#forangsd=$(ls *_dedup_reads.bam)
#for i in $forangsd
#do
#echo $i
#done > samples4angsd.txt

#samtools faidx $TRANSCRIPTOME

### global angsd!
$ANGSD -b samples4angsd.txt -anc $TRANSCRIPTOME -out angsd_global -P 24 -SNP_pval 1e-6 -minMapQ 20 -minQ 20 -setMinDepth 50 -setMaxDepth 6000 -minInd 6 -minMaf 0.01 -GL 1 -doMaf 1 -doMajorMinor 1 -doGlf 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -dosaf 1

### Call angsd on each froggle
# Make a SNP list with major
gunzip angsd_global.mafs.gz
cut -f 1,2,3,4 angsd_global.mafs | tail -n +2 > angsd_global_snplist.txt

# index the sites
$ANGSD sites index angsd_global_snplist.txt
