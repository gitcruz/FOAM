#!/bin/bash

IN=$1;
DB=$2; #/scratch/devel/talioto/ONT/control/dna_cs_fix.fasta
THREADS=$3;
echo "filtering out dna_cs reads from $IN"
b=`basename $IN .fastq`;
module purge;
module load minimap2/2.11
module load gcc/4.9.3-gold xz samtools/1.3
minimap2 -t $THREADS -ax map-ont $DB $IN | samtools calmd -S - $DB > $IN.MINIMAP2.sam
sam_id_multifilter.pl <  $IN.MINIMAP2.sam > $b.clean.fastq 2> $b.contam.fastq
