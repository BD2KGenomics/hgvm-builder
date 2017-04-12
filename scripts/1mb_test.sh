#!/usr/bin/env bash
# 1mb_test.sh: run a 1-megabase test of the vg caller for SVs, manually

set -ex

# For building the graph, what FASTA and VCF should we use?
# VCF also contains SV truth calls
GRCH38_FASTA="../GRCh38.fa"
CHR22_VCF="/cluster/home/charles/SV_HGVM_research/eighth_draft_GRCh38_bks_only_polA25_chrs/ALL.chr22.wgs.integrated_sv_map_v2.20130502.svs.GRCh38.bks.only.genotypes.vcf.gz"

# How close to you have to be to the SV you are looking for, in start coordinates, to count?
THRESHOLD=100

# What region should we process?
REGION_START=36000000
REGION_END=37000000

mkdir -p 1mb

if [ ! -e 1mb/1mb.vg ]; then

    # Make a graph of a 1 megabase region, which contains some SVs in NA12878
    time vg construct -r ${GRCH38_FASTA} -v ${CHR22_VCF} -f -n 22=CM000684.2 -R "22:${REGION_START}-${REGION_END}"  > 1mb/1mb.vg
    
fi

if [ ! -e 1mb/1mb.xg ] || [ ! -e 1mb/1mb.gcsa ]; then

    # Index the graph
    /usr/bin/time -v vg index -x 1mb/1mb.xg -g 1mb/1mb.gcsa 1mb/1mb.vg -k 16
    
fi

if [ ! -e 1mb/NA12878.part.R1.fastq ] || [ ! -e 1mb/NA12878.part.R2.fastq ]; then
    # Go get the relevant reads
    samtools view -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/illumina_platinum_pedigree/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram "chr22:${REGION_START}-${REGION_END}" > 1mb/NA12878.part.bam
    samtools sort -n 1mb/NA12878.part.bam -o 1mb/NA12878.part.byname.bam
    samtools view 1mb/NA12878.part.byname.bam | ./scripts/smartSam2Fastq.py --fq1 1mb/NA12878.part.R1.fastq --fq2 1mb/NA12878.part.R2.fastq --drop_secondary
fi

if [ ! -e 1mb/NA12878.gam ]; then

    # Do the mapping
    time vg map -x 1mb/1mb.xg -g 1mb/1mb.gcsa -f 1mb/NA12878.part.R1.fastq -f 1mb/NA12878.part.R2.fastq > 1mb/NA12878.gam

fi

if [ ! -e 1mb/NA12878.vgpu ]; then

    # Do the pileup
    time vg pileup 1mb/1mb.vg 1mb/NA12878.gam > 1mb/NA12878.vgpu
    
fi

# Make variant calls
time vg call -v -p 1mb/1mb.vg 1mb/NA12878.vgpu -r CM000684.2 -o ${REGION_START} --aug-graph 1mb/aug.vg > 1mb/NA12878.vcf
vg index -x 1mb/aug.xg 1mb/aug.vg

# Grab what we expect
bcftools view -s NA12878 ${CHR22_VCF} -r "22:${REGION_START}-${REGION_END}" | bcftools filter -e 'GT=="0/0" || GT=="0|0" || GT=="0" || FILTER!="PASS"' | grep -v "#" | cut -f1-3,10 > 1mb/want.tsv

# Grab the SVs we actually have. Omit IDs since they're long and useless
cat 1mb/NA12878.vcf | bcftools filter -e 'GT=="0/0" || GT=="0|0" || GT=="0" || FILTER!="PASS" || (SVLEN>-25 && SVLEN<25)' | grep -v "#" | sort -n -k2 | cut -f1-2,10 > 1mb/got.tsv

# Now find out which wanted SVs have a found SV sufficiently close
# See <http://stackoverflow.com/a/1521603>
true >1mb/found.tsv
exec 4<1mb/want.tsv
while read -u4 LINE ; do
    # Where is this SV?
    POSITION=$(echo "${LINE}" | cut -f2)
    # Find the closest value in column 2 of the found SVs
    # See <http://stackoverflow.com/a/17853270>
    NEAREST=$(awk -v c=2 -v t=${POSITION} 'NR==1{d=$c-t;d=d<0?-d:d;v=$c;next}{m=$c-t;m=m<0?-m:m}m<d{d=m;v=$c}END{print v}' 1mb/got.tsv)
    
    if [ "$((POSITION-NEAREST))" -lt "${THRESHOLD}" -a "$((NEAREST-POSITION))" -lt "${THRESHOLD}" ]; then
        # It's in range
        echo "${LINE} matched at ${NEAREST}" >>1mb/found.tsv
    fi
done

# Count up score
EXPECTED=$(cat 1mb/want.tsv | wc -l)
FOUND=$(cat 1mb/found.tsv | wc -l)

PORTION=$(echo "${FOUND} / ${EXPECTED}" | bc -l)

echo "Portion of SVs recalled: ${PORTION}"













