# Population genetics analyses

This is our first-pass analysis looking at population genetics between coastal and inland frogs. We are mostly curious about genomic divergence between these two populations. If we see evidence of divergence and population structure, we will conduct more fine scale analyses.


First, a script that will trim reads, align trimmed, reads to the reference transcriptome, produce sorted bam files for each sample, and finally call sites. Things of note:

1. It only uses the reads in which forward and reverse reads survive.
2. It calls every site, not just variant sites.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=40
#SBATCH --mem 300Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH --output variantcalling.log

# variables
ASSEMBLY=$"Hyla_cinerea_ORP.2.2.8.ORP.fasta"
SAMPLES=$(find data -type f -name "*R1_001.fastq.gz" -printf '%f\n' | sed "s/_R1_001.fastq.gz//g")

### setup
mkdir intermediate_files/
mkdir intermediate_files/trimmed_reads
mkdir intermediate_files/bams/
mkdir intermediate_files/vcfs/

# begin
echo indexing assembly with bwa
bwa index $ASSEMBLY

echo indexing with samtools
samtools faidx $ASSEMBLY

echo starting pipeline
for sample in $SAMPLES
do

echo trimming $sample
trimmomatic PE -threads 40 -baseout intermediate_files/trimmed_reads/$sample.fq.gz \
data/${sample}_R1_001.fastq.gz data/${sample}_R2_001.fastq.gz \
LEADING:3 TRAILING:3 ILLUMINACLIP:barcodes.fa:2:30:10 MINLEN:25

echo aligning $sample
##### NOTE: this only uses reads that are paired, ignores unpaired reads.

bwa mem -t 40 $ASSEMBLY intermediate_files/trimmed_reads/${sample}_1P.fq.gz intermediate_files/trimmed_reads/${sample}_2P.fq.gz \
| samtools view -@20 -Sb - \
| samtools sort -T "$sample" -O bam -@20 -l9 -m2G -o intermediate_files/bams/"$sample".sorted.bam -

echo indexing $sample
samtools index intermediate_files/bams/"$sample".sorted.bam

echo finished with sample prep for $sample
done

echo preparing list of bam files to convert to vcf
find intermediate_files/bams -type f -name "*bam" | tee samplebammies.txt

echo sanity check, are there the correct number of files
if [ $(wc -l samplebammies.txt | cut -d " " -f 1) -eq 31 ]
then
      echo There are the correct number of sample bam files
else
      echo You fucked up somewhere and there are not the correct number of sample bam files
fi


echo converting to vcf format
bcftools mpileup -Ou --threads 40 -f $ASSEMBLY --min-MQ 30 --ignore-RG --max-depth 1000 --bam-list samplebammies.txt | bcftools call --threads 40 --variants-only -m -Ov -o intermediate_files/vcfs/H_cinerea.variantsonly.bcf
```


Then, I applied fairly stringent filtering, since the sample size is low and it is RNA seq.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=40
#SBATCH --mem 300Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH --output variantfiltering.log

VCF="intermediate_files/vcfs/H_cinerea.bcf"
OUT="intermediate_files/vcfs/H_cinerea.filtered.vcf.gz"

# description
Description="

We are applying fairly stringent filtering to our RNAseq dataset. Here is the important information: 


minor allele frequency (maf) = 0.1 
missing data (max-missing) = allow 10% missing data 
minimum quality (minQ) = 30 
minimum depth for a genotype, below which an individual is marked as having missing data (minDP) = 10 
minimum average depth for a site, below which is it marked as missing data (--min-meanDP) = 10 


"

printf "%s" "$Description"


## Filter vcf print
vcftools --bcf $VCF \
--maf 0.1 \
--max-missing 0.9 \
--minQ 30 \
--min-meanDP 10  \
--minDP 10 \
--recode --stdout \
| gzip -c > $OUT
```

Next up, what is the genomic divergence (Fst) between coastal and inland frogs?

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=40
#SBATCH --mem 300Gb
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118
#SBATCH --output fstcalc.log


VCF="intermediate_files/vcfs/H_cinerea.filtered.vcf.gz"

bcftools query -l $VCF

#### Produce population files:
# extract sample names for coastal
bcftools query -l $VCF | grep -E "BOD|CSI|DQ" >  coastalsamples.txt
# extract sample names for inland
bcftools query -l $VCF | grep -E "PWL|WHF|LO|BL" > inlandsamples.txt


#### Use these population files to calculate fst:
vcftools --zvcf ${VCF} \
--weir-fst-pop coastalsamples.txt \
--weir-fst-pop inlandsamples.txt \q
--out coastal-inland
```

We obviously have to look at genetic structure. Run ADMIXTURE.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J refadmix
#SBATCH --output admixture.log
#SBATCH --cpus-per-task=24
#SBATCH --exclude node117,node118

## prep
mkdir intermediate_files/ADMIXTURE
cd intermediate_files/ADMIXTURE

VCF="intermediate_files/vcfs/H_cinerea.filtered.vcf.gz"
OUT="H_cinerea"

plink --vcf $VCF --make-bed --out $OUT 

# run admixture; I have 2 populations, so do up to K = 6
for K in 1 2 3 4 5 6
do
admixture -j40 --cv $OUT.bed $K | tee admix_log${K}.out

# create CV value file:
CV=$(grep "CV error" admix_log${K}.out | cut -d : -f 2 | sed "s/ //g")
printf "%s\t%s\n" "$K" "$CV" >> CV_values.tab
done
```

And then make a fancy structure plot:

_There will eventually be R code here_





