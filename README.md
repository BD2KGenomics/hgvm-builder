# hgvm-builder

This repository contains (or will contain) a Toil script that can build a Human Genome Variation Map from a variety of data sources.

The basic workflow will be:

1. Start with a GRC assembly, like ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.24_GRCh38.p9/GCA_000001405.24_GRCh38.p9_assembly_structure/

2. Build it into a HAL-fromat graph with Cactus

3. Convert the graph into vg format

4. Add variation to the graph from one or more VCF files

5. Index the graph for alignment and variant calling

6. Align, call, and integrate collections of short read data

## Data Sources

Data from a variety of sources will be combined to build the HGVM.

### HGVM 0.1 Alpha

For the first pass, only two data sources will be incorporated.

- [ ] The GRCh38 assembly, including alt loci

- [ ] The 1000 Genomes Project small variant VCFs

- [ ] The 1000 Genomes Project structural variant VCFs

### Subsequent Sprints

Once the basic infrastructure for building a graph, the following data sources will be added.

- [ ] GRCh38 novel and fix patches

- [ ] The Simons Genome Diversity Project

- [ ] The Aftican Genome Variation Project

- [ ] The Personal Genome Project

- [ ] Platinum Genomes

- [ ] The Wash-U PacBia Assemblies




