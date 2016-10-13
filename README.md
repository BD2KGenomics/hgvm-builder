# hgvm-builder

This repository contains (or will contain) a Toil script that can build a Human Genome Variation Map from a variety of data sources.

The basic workflow will be:

1. Start with a GRC assembly, like ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.24_GRCh38.p9/GCA_000001405.24_GRCh38.p9_assembly_structure/

2. Break out top-level contigs (chromosomes and unplaced scaffolds) and construct them into graphs, along with applicable VCFs, using vg

3. Merge these graphs together and index them

4. Align and merge in alt locus sequences, fix and novel patches and other full-length assemblies.

5. Align, call, and integrate collections of short read data
