#!/bin/bash

# coastal Hyla cinerea transcriptome assembly and downstream bash work
# start in molly/ folder


### Prepare coastal samples for assembly

cd coastal/ 

# first, concatenate reads by coastal treatment (0-0, 0-4, 4-4, 4-6)
for i in reads/0_0_*R1_001.fastq; do cat $i; echo; done > coastal.zerozero.R1.fastq 
for f in reads/0_0_*R2_001.fastq; do cat $f; echo; done > coastal.zerozero.R2.fastq


for g in reads/0_6_*R1_001.fastq; do cat $g; echo; done > coastal.zerosix.R1.fastq 
for h in reads/0_6_*R2_001.fastq; do cat $h; echo; done > coastal.zerosix.R2.fastq 


for l in reads/4_4_*R1_001.fastq; do cat $l; echo; done > coastal.fourfour.R1.fastq 
for m in reads/4_4_*R2_001.fastq; do cat $m; echo; done > coastal.fourfour.R2.fastq


for n in reads/4_6_*R1_001.fastq; do cat $n; echo; done > coastal.foursix.R1.fastq 
for o in reads/4_6_*R2_001.fastq; do cat $o; echo; done > coastal.foursix.R2.fastq 


# subsample each treatment to equal coverage. since there are 4 treatments, we are sampling each *treatment group* to 12.5M reads, and then running a transcriptome assembly on the concatenated dataset of 50M reads
seqtk sample -s100 coastal.zerozero.R1.fastq 12500000 > subsamp.coastal.zerozero.R1.fastq 
seqtk sample -s100 coastal.zerozero.R2.fastq 12500000 > subsamp.coastal.zerozero.R2.fastq 

seqtk sample -s100 coastal.zerosix.R1.fastq 12500000 > subsamp.coastal.zerosix.R1.fastq 
seqtk sample -s100 coastal.zerosix.R2.fastq 12500000 > subsamp.coastal.zerosix.R2.fastq 

seqtk sample -s100 coastal.fourfour.R1.fastq 12500000 > subsamp.coastal.fourfour.R1.fastq 
seqtk sample -s100 coastal.fourfour.R2.fastq 12500000 > subsamp.coastal.fourfour.R2.fastq 

seqtk sample -s100 coastal.foursix.R1.fastq 12500000 > subsamp.coastal.foursix.R1.fastq 
seqtk sample -s100 coastal.foursix.R2.fastq 12500000 > subsamp.coastal.foursix.R2.fastq 


# concatenate forward and reverse reads, each should 50 million reads

for i in subsamp.coastal*R1.fastq; do cat $i; echo; done > subsamp.coastal.R1.fastq 
for f in subsamp.coastal*R2.fastq; do cat $f; echo; done > subsamp.coastal.R2.fastq

# delete intermediate files
rm coastal.zero* coastal.four* subsamp.coastal.zero* subsamp.coastal.four*


### Prepare inland samples for assembly
cd ../inland/

# first, concatenate reads by inland treatment (0-0, 0-4, 4-4, 4-6)
for i in reads/0_0_*R1_001.fastq; do cat $i; echo; done > inland.zerozero.R1.fastq 
for f in reads/0_0_*R2_001.fastq; do cat $f; echo; done > inland.zerozero.R2.fastq 


for g in reads/0_6_*R1_001.fastq; do cat $g; echo; done > inland.zerosix.R1.fastq 
for h in reads/0_6_*R2_001.fastq; do cat $h; echo; done > inland.zerosix.R2.fastq 


for l in reads/4_4_*R1_001.fastq; do cat $l; echo; done > inland.fourfour.R1.fastq 
for m in reads/4_4_*R2_001.fastq; do cat $m; echo; done > inland.fourfour.R2.fastq 


for n in reads/4_6_*R1_001.fastq; do cat $n; echo; done > inland.foursix.R1.fastq 
for o in reads/4_6_*R2_001.fastq; do cat $o; echo; done > inland.foursix.R2.fastq 


# subsample each treatment to equal coverage. since there are 4 treatments, we are sampling each *treatment group* to 12.5M reads, and then running a transcriptome assembly on the concatenated dataset of 50M reads
seqtk sample -s100 inland.zerozero.R1.fastq 12500000 > subsamp.inland.zerozero.R1.fastq 
seqtk sample -s100 inland.zerozero.R2.fastq 12500000 > subsamp.inland.zerozero.R2.fastq 

seqtk sample -s100 inland.zerosix.R1.fastq 12500000 > subsamp.inland.zerosix.R1.fastq 
seqtk sample -s100 inland.zerosix.R2.fastq 12500000 > subsamp.inland.zerosix.R2.fastq 

seqtk sample -s100 inland.fourfour.R1.fastq 12500000 > subsamp.inland.fourfour.R1.fastq 
seqtk sample -s100 inland.fourfour.R2.fastq 12500000 > subsamp.inland.fourfour.R2.fastq 

seqtk sample -s100 inland.foursix.R1.fastq 12500000 > subsamp.inland.foursix.R1.fastq 
seqtk sample -s100 inland.foursix.R2.fastq 12500000 > subsamp.inland.foursixe.R2.fastq 

# concatenate forward and reverse reads, each should be 50 million reads

for i in subsamp.inland*R1.fastq; do cat $i; echo; done > subsamp.inland.R1.fastq 
for f in subsamp.inland*R2.fastq; do cat $f; echo; done > subsamp.inland.R2.fastq

# delete intermediate files
rm inland.zero* inland.four* subsamp.inland.zero* subsamp.inland.four*


### Combine reads from inland and coastal treatments

# merge subsampled datasets into a single file
cat */*R1.fastq > combined.subsamp.R1.fastq
cat */*R2.fastq > combined.subsamp.R2.fastq


### Run the Oyster River Protocol to assemble a transcriptome
/home/summersk/programs/Oyster_River_Protocol/oyster35_55.mk main \
MEM=750 \
CPU=30 \
READ1=combined.subsamp.R1.fastq \
READ2=combined.subsamp.R2.fastq \
RUNOUT=100Mrun



cp assemblies/100Mrun.orthomerged.fasta .

# Rename all the transcripts because aesthetics.
awk '/^>/{print ">Transcript_" ++i; next}{print}'  100Mrun.orthomerged.fasta > 100Mrun.cinerea.fasta

###### Transcriptome assembly complete


### Trim and clean sample reads

# clean coastal first
cd coastal/

# Remove adaptors and trim the reads from each sample, first get the samples
samples=$(ls reads/*fastq | sed "s/R1_001.fastq//g" | sed "s/R2_001.fastq//g" | uniq)

# Run trimmomatic
cd reads/
(ls *fastq | sed "s/_R1_001.fastq//g" | sed "s/_R2_001.fastq//g" | uniq) | \
parallel -j 5 trimmomatic-0.36.jar PE -threads 4 \
-baseout /home/summersk/molly/coastal/rcorr/{}.TRIM.fastq {}_R1_001.fastq {}_R2_001.fastq \
LEADING:3 TRAILING:3 ILLUMINACLIP:barcodes.fa:2:30:10 MINLEN:25 

# R corrector  
cd ..
(ls rcorr/*P.fastq | sed "s/.TRIM_1P.fastq//g" | sed "s/.TRIM_2P.fastq//g"  | uniq | grep -v coastalsubsamp) | \
parallel -j 6 run_rcorrector.pl -t 4 -k 31 -1 {}.TRIM_1P.fastq -2 {}.TRIM_2P.fastq 

mv *cor.fq rcorr/



### clean inland second
cd ../inland/

# Remove adaptors and trim the reads from each sample, first get the samples
samples=$(ls reads/*fastq | sed "s/R1_001.fastq//g" | sed "s/R2_001.fastq//g" | uniq)

# Run trimmomatic
cd reads
(ls *fastq | sed "s/_R1_001.fastq//g" | sed "s/_R2_001.fastq//g" | uniq) | \
parallel -j 5 trimmomatic-0.36.jar PE -threads 4 \
-baseout /home/summersk/molly/inland/rcorr/{}.TRIM.fastq {}_R1_001.fastq {}_R2_001.fastq \
LEADING:3 TRAILING:3 ILLUMINACLIP:barcodes.fa:2:30:10 MINLEN:25 

# R corrector  
cd ..
(ls rcorr/*P.fastq | sed "s/.TRIM_1P.fastq//g" | sed "s/.TRIM_2P.fastq//g"  | uniq | grep -v inlandsubsamp) | \
parallel -j 6 run_rcorrector.pl -t 4 -k 31 -1 {}.TRIM_1P.fastq -2 {}.TRIM_2P.fastq 

mv *cor.fq rcorr/

cd ..


### Read mapping

# Pseudo-quantification with kallisto, build the index first
kallisto index -i cinerea.idx  100Mrun.cinerea.fasta

# This will list all the samples
cd inland/rcorr/
inlsamples=$(ls *P.cor.fq | sed "s/.TRIM_1P.cor.fq//g" | sed "s/.TRIM_2P.cor.fq//g" | uniq  | grep -v subsamp)

cd ../../coastal/rcorr
coastsamples=$(ls *P.cor.fq | sed "s/.TRIM_1P.cor.fq//g" | sed "s/.TRIM_2P.cor.fq//g" | uniq  |
grep -v subsamp)

cd ../../

# Make directories for each sample:
mkdir kallisto_quants
cd kallisto_quants
for i in $inlsamples; do mkdir $i; done
for i in $coastsamples; do mkdir $i; done
cd ..


# Perform the actual pseudo-quantification for all of the samples + technical replicates

parallel -j 10 kallisto quant -i /home/summersk/molly/cinerea.idx -o kallisto_quants/{} -b 100 \
inland/rcorr/{}.TRIM_1P.cor.fq inland/rcorr/{}.TRIM_2P.cor.fq ::: $inlsamples

parallel -j 16 kallisto quant -i /home/summersk/molly/cinerea.idx -o kallisto_quants/{} -b 100 \
coastal/rcorr/{}.TRIM_1P.cor.fq coastal/rcorr/{}.TRIM_2P.cor.fq ::: $coastsamples


#### Annotation with diamond

diamond blastx -d /home/summersk/peptide_databases/xen.dmnd -q 100Mrun.cinerea.fasta -o cinerea2xen.m8 --threads 32 
# 31,442 queries aligned

# sort by top hit
sort cinerea2xen.m8 -k 1,1 -k11,11g | sort -u -k 1,1 --merge > cinerea_Xen95_annotation.tsv


