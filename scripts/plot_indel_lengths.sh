#!/usr/bin/env bash
# plot_indel_lengths.sh: Make an indel length plot for predicted and true VCFs

set -ex

BUILD_DIR="builds/2017-04-25"
TRUTH_VCF="/cluster/home/charles/SV_HGVM_research/eighth_draft_GRCh38_bks_only_polA25_chrs/ALL.chr22.wgs.integrated_sv_map_v2.20130502.svs.GRCh38.bks.only.genotypes.vcf.gz"
CONTIG="22"
CALLED_VCF="builds/2017-04-24/eval/sv/NA12878/calls.vcf"
#CALLED_VCF="builds/2017-04-20/chr22/recall/calls.vcf"
# How close to you have to be to the SV you are looking for, in start coordinates, to count?
THRESHOLD=100
SCRATCH_DIR="scratch"
mkdir -p "${SCRATCH_DIR}"

true > "${SCRATCH_DIR}/svlens.tsv"
true > "${SCRATCH_DIR}/expected.tsv"

for SAMPLE in NA12878 NA12889 NA12890; do

    # Find the calls for this sample
    CALL_VCF="${BUILD_DIR}/eval/sv/${SAMPLE}/calls.vcf"
    cat "${CALL_VCF}" | bcftools filter -e 'GT=="0/0" || GT=="0|0" || GT=="0" || FILTER!="PASS" || (SVLEN>-25 && SVLEN<25)' | grep -v "^#" | cut -f8 | cut -d';' -f2 | cut -d'=' -f2 | tr ',' '\n' >> "${SCRATCH_DIR}/svlens.tsv"

    # And the truth
    bcftools view -s "${SAMPLE}" "${TRUTH_VCF}" -r "${CONTIG}" | bcftools filter -e 'GT=="0/0" || GT=="0|0" || GT=="0"' | ./scripts/indelLengths.py --distinguish >> "${SCRATCH_DIR}/expected.tsv"
    
done

    ./scripts/histogram.py "${SCRATCH_DIR}/svlens.tsv" "${SCRATCH_DIR}/expected.tsv" \
        --title "Deletion lengths" \
        --x_label "Length (bp)" --y_label "Deletion count" --save svlens.del.png \
        --line --no_n \
        --bins 30 --no_zero_ends --split_at_zero --x_max -25 \
        --legend_overlay "upper left" \
        --style "-" \
        --category_labels "Detected" "Expected" \
        --width 6 --height 4

    ./scripts/histogram.py "${SCRATCH_DIR}/svlens.tsv" "${SCRATCH_DIR}/expected.tsv" \
        --title "Insertion lengths" \
        --x_label "Length (bp)" --y_label "Insertion count" --save svlens.ins.png \
        --line --no_n \
        --bins 30 --no_zero_ends --split_at_zero --x_min 25 \
        --legend_overlay "upper right" \
        --style "-" \
        --category_labels "Detected" "Expected" \
        --width 6 --height 4
        
