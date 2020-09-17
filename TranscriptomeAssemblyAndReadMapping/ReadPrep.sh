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
cat */*R1.fastq > subsamp.R1.fastq
cat */*R2.fastq > subsamp.R2.fastq
