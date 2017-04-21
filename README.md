# hgvm-builder

This repository contains a Toil script that can build a Human Genome Variation Map from a variety of data sources, and evaluate it.

The basic workflow is:

1. Start with a vg graph already built, or a HAL file built from an assembly like [GRCh38](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.24_GRCh38.p9/GCA_000001405.24_GRCh38.p9_assembly_structure/)

2. Convert into a merged vg graph

3. Add variation to the graph from one or more VCF files

4. Index the graph for alignment and variant calling

5. Align smaples and call variants

6. Evaluate calls with truth VCFs or the realignment of assembly fragments

## Installation

Install (for development) with a `python setup.py develop`. This should also install Python dependencies.

Run a tiny example with `scripts/test_builder.sh` to make sure everything is working.

Then use the `build-hgvm` tool.

## Azure Deployment

To run on Azure, first make sure you have your input data in Azure Storage, and that you can connect to it:

```
export AZURE_STORAGE_CONNECTION_STRING='DefaultEndpointsProtocol=https;AccountName=YOURNAMEHERE;AccountKey=YOURKEYHERE'
```

Use the `scripts/toilAzureCleanup.sh` script to clear out old undeleted job stores that Toil may leave lying about.

Upload your data to Azure with copy-hgvm:

```
PYTHONPATH="src" python -m hgvmbuilder.parallelcopy ./tree file:./things azure:account:container/path/to/things
```

## Project Roadmap

Data from a variety of sources will be combined to build the HGVM.

### HGVM 0.1 Alpha

For the first pass, only a few data sources will be incorporated.

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




