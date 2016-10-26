# hgvm-builder

This repository contains (or will contain) a Toil script that can build a Human Genome Variation Map from a variety of data sources.

The basic workflow will be:

1. Start with a GRC assembly, like ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.24_GRCh38.p9/GCA_000001405.24_GRCh38.p9_assembly_structure/

2. Break out top-level contigs (chromosomes and unplaced scaffolds) and construct them into graphs, along with applicable VCFs, using vg

3. Merge these graphs together and index them

4. Align and merge in alt locus sequences, fix and novel patches and other full-length assemblies.

5. Align, call, and integrate collections of short read data

## Data Sources

Data from a variety of sources will be combined to build the HGVM.

### Sprint 1

For the first pass, only two data sources will be incorporated.

- [ ] The GRCh38 assembly, including alt loci and novel and fix patches.

- [ ] The 1000 Genomes Project lifted-over VCFs

### Subsequent Sprints

Once the basic infrastructure for building a graph, the following data sources will be added.

- [ ] The Simons Genome Diversity Project

- [ ] The Aftican Genome Variation Project

- [ ] The Personal Genome Project

- [ ] Platinum Genomes

- [ ] The Wash-U PacBia Assemblies




