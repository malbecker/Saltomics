# Population genetics analyses

This is our first-pass analysis looking at population genetics between coastal and inland frogs. We are mostly curious about genomic divergence between these two populations. If we see evidence of divergence and population structure, we will conduct more fine scale analyses.


First, a script that will trim reads, align trimmed, reads to the reference transcriptome, produce sorted bam files for each sample, and finally call sites. Of note: This script only uses the reads in which forward and reverse reads survive.

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

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node117,node118
#SBATCH --output variantfiltering.log

VCF="intermediate_files/vcfs/H_cinerea.variantsonly.vcf"
OUT="intermediate_files/vcfs/H_cinerea.filtered.vcf"

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


## Filter vcf
vcftools --vcf $VCF \
--maf 0.1 \
--max-missing 0.9 \
--minQ 20 \
--remove-indels \
--recode --stdout > $OUT

# --min-meanDP 10  \
# --minDP 10 \
```

Next up, what is the genomic divergence (Fst) between coastal and inland frogs? Even though sample sizes are super low for individual comparisons, we make those as well as a sanity check.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH --ntasks=24
#SBATCH --exclude=node117,node118
#SBATCH --output fstcalc.log


VCF="intermediate_files/vcfs/H_cinerea.filtered.vcf"  # change to filtered file later.
OUTDIR="intermediate_files/fsts"

mkdir intermediate_files/fsts

##########################################
######### while read all  ################
##########################################

COMPARISONS="FstComparisons.txt"

sed 1d $COMPARISONS | while IFS=$'\t' read -r comparison pop1 searchterm1 pop2 searchterm2
do
        echo $comparion with $searchterm1 and $searchterm2
        echo sample prep
        bcftools query -l $VCF | grep $pop1 > ${OUTDIR}/pop1samples.txt
        bcftools query -l $VCF | grep $pop2 > ${OUTDIR}/pop2samples.txt
        echo run fst calculations for $comparison
        vcftools --vcf $VCF --weir-fst-pop ${OUTDIR}/pop1samples.txt --weir-fst-pop ${OUTDIR}/pop2samples.txt --out ${OUTDIR}/$comparison > ${OUTDIR}/$comparison.fst.log 2>&1
done


##### coastal inland comparison

bcftools query -l $VCF | grep -E "BOD|CSI|DQ" > ${OUTDIR}/pop1samples.txt
bcftools query -l $VCF | grep -E "PWL|WHF|LO|BL" > ${OUTDIR}/pop2samples.txt
echo run fst calculations for coastal vs inland
vcftools --vcf $VCF --weir-fst-pop ${OUTDIR}/pop1samples.txt --weir-fst-pop ${OUTDIR}/pop2samples.txt --out ${OUTDIR}/coastal-inland > ${OUTDIR}/coastal-inland.fst.log 2>&1

### extract fst values:

printf "Comparison\tPopulation_1\tPopulation_2\tMean_Fst\tWeighted_Fst\tSampleSize\n" > Fstvalues.tsv

comps=$(ls ${OUTDIR}/*fst.log)

for comp in $comps
do
pop1=$(grep "\\-\\-out" $comp | cut -f 2 -d " " | cut -d- -f1 | sed "s/^.*\\///g")
pop2=$(grep "\\-\\-out" $comp | cut -f 2 -d " " | cut -d- -f2)
fst=$(grep "mean Fst" $comp | cut -d: -f2 | sed "s/ //g")
weir=$(grep "weighted" $comp | cut -d: -f2 | sed "s/ //g")
N=$(grep "Individuals" $comp | cut -f 4 -d " ")
comppops=$(echo $comp | sed "s/^.*\\///g" | sed "s/.fst.log//g")
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$comppops" "$pop1" "$pop2" "$fst" "$weir" "$N" >> Fstvalues.tsv
done

## cleanup
rm ${OUTDIR}/pop1samples.txt
rm ${OUTDIR}/pop2samples.txt
rm ${OUTDIR}/*fst.log

```

We obviously have to look at genetic structure. We want to run ADMIXTURE, but first we have to convert our data to a usable format via `plink2`. Because plink only accepts a certain number of contigs and we are over that, we have to do a bit of nonsense.

```bash
cd intermediate_files/vcfs/

grep "^##contig=<ID" H_cinerea.filtered.vcf | cut -f3 -d= | cut -f1 -d, | sort | uniq > contigs.txt

mkdir tmp_vcfs/
mkdir tmp_plink/

# pull out each contig into its own vcf
for tig in `cat contigs.txt`; do bcftools view H_cinerea.filtered.vcf.gz ${tig} > tmp_vcfs/${tig}.vcf; done

# make each contig vcf into the appropriate bed/bim/fam file format
for tig in `cat contigs.txt`; do /mnt/lustre/macmaneslab/ams1236/software/plink2 --double-id --allow-extra-chr --vcf tmp_vcfs/${tig}.vcf --make-bed --out tmp_plink/${tig}; done

# cat them all together (and hope it works).
cd tmp_plink
cat *bed > tmp
mv tmp > H_cinerea.filtered.bed

cat *bim > tmp
mv tmp > H_cinerea.filtered.bim

cat *fam > tmp
mv tmp > H_cinerea.filtered.fam
```

Now we can run ADMIXTURE.

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





