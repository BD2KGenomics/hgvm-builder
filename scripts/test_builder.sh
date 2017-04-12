#!/usr/bin/env bash

set -e

rm -Rf test_tree test_build
build-hgvm ./test_tree test_build \
    --base_vg_url file:`pwd`/test/tiny-flat.vg \
    --vcf_contig "x" \
    --vcf_url file:`pwd`/test/tiny.vcf \
    --sv_sample_name "1" \
    --sv_sample_fastq_url file:`pwd`/test/tiny.fastq \
    --logInfo \
    --realTimeLogging \
    --sample_fastq_url file:`pwd`/test/tiny.fastq \
    --eval_sequences_url file:`pwd`/test/tiny.seqs \
    --control_graph_url file:`pwd`/test/tiny-flat.vg \
    --control_graph_xg_url file:`pwd`/test/tiny-flat.xg \
    --control_graph_gcsa_url file:`pwd`/test/tiny-flat.gcsa \
