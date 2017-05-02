#!/usr/bin/env bash
# validate_sv_caller.sh: test the vg sv caller to make sure it works
set -e

# TODO: keep files in a temp directory

# Make a graph
vg construct -r test/x.fa > x-flat.vg
vg index -x x-flat.xg -g x-flat.gcsa -k 16 x-flat.vg

# Make a file of characteristic indel sequences, with no SNPs
true > indel_seqs.txt
for SEQ_NUM in {1..10}; do
    # Create a seeded sequence
    SEQ=$(vg sim -n 1 -l 100 -x x-flat.xg -e 0.00 -i 0 -f -s "${SEQ_NUM}")
    if [ "${SEQ_NUM}" -lt 5 ]; then
        # Cut out the middle to make it a deletion
        echo "${SEQ:0:30}${SEQ:70:30}" >> indel_seqs.txt
    else
        # Insert something to make it an insertion
        INSERT_SEQ=$(vg sim -n 1 -l 50 -x x-flat.xg -e 0.00 -i 0 -f -s "${SEQ_NUM}00")
        echo "${SEQ:0:50}${INSERT_SEQ}${SEQ:50:50}" >> indel_seqs.txt
    fi
done

# Map them, making sure to stamp down the ends hard
vg map -x x-flat.xg -g x-flat.gcsa --reads indel_seqs.txt --full-l-bonus 100 -t 1 > svs.gam

# Make a new reference graph
vg mod -i svs.gam x-flat.vg | vg mod --unchop - > reference.vg
vg index -x reference.xg -g reference.gcsa -k 16 reference.vg

# Find its ends
LEFT_ID=$(vg find -x reference.xg -p x:1 | vg view -j - 2>/dev/null | jq '.node[].id')
RIGHT_ID=$(vg find -x reference.xg -p x:1000 | vg view -j - 2>/dev/null | jq '.node[].id')

vg view -dp reference.vg | dot -Tsvg -o test.svg

# Simulate a haplotype pair that touches both end nodes
vg sim -n 1000 -l 1000 -x reference.xg -s 100 -e 0 -i 0 -J | grep "node_id.: ${LEFT_ID}}" | grep "node_id.: ${RIGHT_ID}}" | head -n2 | vg view -JGa - > haplotype.gam

# Plot it (we have some nice overlaps and nesting)
vg view -dp reference.vg -A haplotype.gam | dot -Tsvg -o test.svg

# Decide what variants it has, duplicating it out to ensure calls
cat haplotype.gam haplotype.gam haplotype.gam haplotype.gam haplotype.gam | vg pileup reference.vg - > haplotype.vgpu

# And make it a graph
vg view -aj haplotype.gam | jq -r '">" + .name, .sequence' > haplotype.fa
vg construct -r haplotype.fa > haplotype.vg
vg index -x haplotype.xg haplotype.vg

for READ_LENGTH in 50 40 30 20 10; do
    # Go down in read length until there's a difference

    # Simulate some sample reads from it
    # 100 * 60 / 1000 = 6x coverage ish
    # 1000 * 60 / 1000 = 60x coverage ish
    DESIRED_COVERAGE=60
    REF_LENGTH=1000

    READ_COUNT=$(($DESIRED_COVERAGE * $REF_LENGTH / $READ_LENGTH))
    echo "Making ${READ_COUNT} reads..."
    vg sim -n ${READ_COUNT} -l ${READ_LENGTH} -x haplotype.xg -e 0.00 -i 0.00 -s 9999 > reads.txt

    # Map them
    vg map -x reference.xg -g reference.gcsa --reads reads.txt > reads.gam
    vg pileup reference.vg reads.gam > reads.vgpu

    #vg view -d reference.vg -A reads.gam | dot -Tsvg -o test.svg

    vg call reference.vg haplotype.vgpu | grep -v "0/0" | grep -v "#" | sort -n -k2 | cut -f1,2,10 | cut -f1 -d':' | sed s_1/0_0/1_g > truth.tsv
    vg call reference.vg reads.vgpu | grep -v "0/0" | grep -v "#" | sort -n -k2 | cut -f1,2,10 | cut -f1 -d':' | sed s_1/0_0/1_g > calls.tsv

    echo "Differences at read length ${READ_LENGTH}:"
    # Get diff for calls based off truth.
    # Diff returns 1 if there's a difference.
    diff calls.tsv truth.tsv
done


