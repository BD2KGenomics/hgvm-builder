#!/usr/bin/env bash
# 1mb_test.sh: run a 1-megabase test of the vg caller for SVs, manually

set -e

GRCH38_FASTA="../GRCh38.fa"
CHR22_VCF="/cluster/home/charles/SV_HGVM_research/eighth_draft_GRCh38_bks_only_polA25_chrs/ALL.chr22.wgs.integrated_sv_map_v2.20130502.svs.GRCh38.bks.only.genotypes.vcf.gz"

mkdir -p 1mb

# Make a graph of a 1 megabase region, which contains some SVs in NA12878
time vg construct -r ${GRCH38_FASTA} -v ${CHR22_VCF} -f -n 22=CM000684.2 -R 22:17000000-18000000  > 1mb/1mb.vg

# Index the graph
/usr/bin/time -v vg index -x 1mb/1mb.xg -g 1mb/1mb.gcsa 1mb/1mb.vg -k 16

# Do the mapping
time vg map -x 1mb/1mb.xg -g 1mb/1mb.gcsa -f NA12878.part.R1.fastq -f NA12878.part.R2.fastq > 1mb/NA12878.gam

# Do the pileup
time vg pileup 1mb/1mb.vg 1mb/NA12878.gam > 1mb/NA12878.vgpu

# Make variant calls
time vg call -v -p 1mb/1mb.vg 1mb/NA12878.vgpu -r CM000684.2 > 1mb/NA12878.vcf

# Grab what we expect
bcftools view -s NA12878 ${CHR22_VCF} | cut -f1-3,10 | grep -v '0/0' | grep -v '0|0' | awk  -F $'\t' 'BEGIN {OFS = FS} {if ($2 > 17000000 && $2 < 18000000) { $2 = $2 - 17000000; print $0 }}'
