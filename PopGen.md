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



"

printf "%s" "$Description"


## Filter vcf
vcftools --vcf $VCF \
--maf 0.1 \
--max-missing 0.9 \
--minQ 20 \
--remove-indels \
--recode --stdout > $OUT

```

Next up, what is the genomic divergence (Fst) between coastal and inland frogs? Even though sample sizes are super low for individual comparisons, we make those as well as a sanity check. Here is a table of all the comparisons. This is the file `FstComparisons.txt` referenced in the script below it.

```
Comparison      Pop1    SearchTerm1     Pop2    SearchTerm2
coastal-inland  coastal BOD1|BOD2|CSI|DQ        inland  PWL|WHF|LO|BL
BL-BOD1 BL      BL      BOD1    BOD1
BL-BOD2 BL      BL      BOD2    BOD2
BL-CSI  BL      BL      CSI     CSI
BL-DQ   BL      BL      DQ      DQ
BL-LO   BL      BL      LO      LO
BL-PWL  BL      BL      PWL     PWL
BL-WHF  BL      BL      WHF     WHF
BOD1-BOD2       BOD1    BOD1    BOD2    BOD2
BOD1-CSI        BOD1    BOD1    CSI     CSI
BOD1-DQ BOD1    BOD1    DQ      DQ
BOD1-LO BOD1    BOD1    LO      LO
BOD1-PWL        BOD1    BOD1    PWL     PWL
BOD1-WHF        BOD1    BOD1    WHF     WHF
BOD2-CSI        BOD2    BOD2    CSI     CSI
BOD2-DQ BOD2    BOD2    DQ      DQ
BOD2-LO BOD2    BOD2    LO      LO
BOD2-PWL        BOD2    BOD2    PWL     PWL
BOD2-WHF        BOD2    BOD2    WHF     WHF
CSI-DQ  CSI     CSI     DQ      DQ
CSI-LO  CSI     CSI     LO      LO
CSI-PWL CSI     CSI     PWL     PWL
CSI-WHF CSI     CSI     WHF     WHF
DQ-LO   DQ      DQ      LO      LO
DQ-PWL  DQ      DQ      PWL     PWL
DQ-WHF  DQ      DQ      WHF     WHF
LO-PWL  LO      LO      PWL     PWL
LO-WHF  LO      LO      WHF     WHF
PWL-WHF PWL     PWL     WHF     WHF
```

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

# Now we have to concatenate them all together.
mkdir plink_merged/
# just the first fam file:
first=$(head -n1 contigs.txt)
cat tmp_plink/${first}.fam > plink_merged/H_cinerea.filtered.merged.fam

# for loop for the bim files:
for tig in `cat contigs.txt`; do cat tmp_plink/${tig}.bim; done > plink_merged/H_cinerea.filtered.merged.bim

# for loop for the bed files:
(echo -en "\x6C\x1B\x01"; for tig in `cat contigs.txt`; do tail -c +4 tmp_plink/${tig}.bed; done) > plink_merged/H_cinerea.filtered.merged.bed

# apparently ADMIXTURE only works with standard chromosome names (integers, X, Y, etc). So here is some more bullshit that changes column 1 to 0 (unknown chromosome), then changes column 2 to what was column 1 (ie the contig id). This then solves our issues with too many contig IDs, etc and we can get ADMIXTURE to run, finally.
awk {'printf ("0\t%s\t%s\t%s\t%s\t%s\t\n", $1, $3, $4, $5, $6)'} H_cinerea.filtered.merged.bim > tmp
mv tmp H_cinerea.filtered.merged.bim 

```

Now we can run ADMIXTURE.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J admix
#SBATCH --output admixture.log
#SBATCH --cpus-per-task=24
#SBATCH --exclude node117,node118

## prep
mkdir intermediate_files/ADMIXTURE

BED="intermediate_files/vcfs/plink_merged/H_cinerea.filtered.merged.bed"

# admixture is in:
#module load linuxbrew/colsa
# run admixture; I have 2 populations, so do up to K = 6
for K in 1 2 3 4 5 6
do
admixture -j24 --cv $BED $K | tee intermediate_files/ADMIXTURE/admix_log${K}.out

# create CV value file:
CV=$(grep "CV error" intermediate_files/ADMIXTURE/admix_log${K}.out | cut -d : -f 2 | sed "s/ //g")
printf "%s\t%s\n" "$K" "$CV" >> intermediate_files/ADMIXTURE/CV_values.tab
done
```

And then make a fancy structure plot:

_There will eventually be R code here_





