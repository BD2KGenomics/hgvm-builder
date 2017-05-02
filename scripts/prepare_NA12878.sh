#!/usr/bin/env bash
# prepare_NA12878.sh: download NA12878 reads and prepare chr22-related reads

# We need to include all the reads on chr22 or an associated alt, but also any
# reads that are pair partners with those reads. We will keep them all in paired
# FASTQs.

set -ex

# Download the CRAM
wget --progress=dot:giga ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/illumina_platinum_pedigree/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram -O NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram

# Sort reads by name, so pair partners are together
samtools sort -n NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram -o NA12878.byname.bam

# Go through in name order and collect pairs where one end touches chr22 or a related alt
samtools view NA12878.byname.bam | awk '{if ($3 ~ /chr22(_.*)?$/ || $7 ~ /chr22(_.*)?$/) print}' | ./scripts/smartSam2Fastq.py --fq1 NA12878.chr22related.R1.fastq --fq2 NA12878.chr22related.R2.fastq --drop_secondary --expect_paired
