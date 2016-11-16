#!/usr/bin/env python
# hgvm-builder build.py: Command-line tool to build a Human Genome Variation Map
"""

Takes a GRC-format assembly strucutre, and a VCF, and produces a vg graph of the
assembly (including novel and fix patches, alternate loci, and
unlocalized/poorly localized scaffonds) including the variation from the VCF.

Parallelizes using the Toil system.

Requires the "vg" binary to be available on the PATH on all nodes, or provided
in a static build usable by all Toil nodes at a URL.

"""

import argparse
import logging
import os
import os.path
import sys
import subprocess
import gzip
import shutil

import toil
from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from .plan import ReferencePlan
from . import grcparser
from . import thousandgenomesparser
from .transparentunzip import TransparentUnzip

# Get a submodule-global logger
Logger = logging.getLogger("build")

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("--assembly_url",
        default=("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCA_000001405.24_GRCh38.p9/"
        "GCA_000001405.24_GRCh38.p9_assembly_structure"),
        help="root of input assembly structure in GRC format")
    parser.add_argument("--vcfs_url", default=("ftp://ftp.1000genomes.ebi.ac.uk/"
        "vol1/ftp/release/20130502/supporting/GRCh38_positions"),
        help="directory of VCFs per chromosome") 
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def create_plan(assembly_url, vcfs_url):
    """
    Given an FTP or file url to the root of a GRC-format assembly_structure
    directory tree, and an FTP or file URL to a directory of chrXXX VCFs,
    produce a ReferencePlan describing that assembly and those VCFs.
    """
    
    # Make the plan
    plan = ReferencePlan()

    # Parse the assembly and populate the plan    
    grcparser.parse(plan, assembly_url)
    
    # Parse the VCF directory and add the VCFs
    thousandgenomesparser.parse(plan, vcfs_url)
    
    # Return the completed plan
    return plan
    
def main_job(job, options, plan):
    """
    Root Toil job. Execute the plan. Returns a list of graph file IDs to export.
    """

    # Make sure we have the Toil IDs for all the top-level chromosome FASTAs.
    fasta_ids = list(plan.for_each_primary_fasta())
    
    # Grab the VCF names and file IDs
    vcf_ids = dict(plan.for_each_vcf_id_by_chromosome())
    
    # Build a list of VG graph promises
    contig_graph_ids = []
    
    RealtimeLogger.info("Building graph from {} primary FASTAs and {} "
        "VCFs".format(len(fasta_ids), len(vcf_ids)))
    
    # Build in parallel for each chromosome
    for chromosome_name, vcf_id in vcf_ids.iteritems():
        # For each chromosome and its VCF
        
        # Make a job to build a graph with that VCF and the relevant FASTA
        
        # TODO: for now we'll just globally replicate the FASTAs and let
        # cacheing sort it out. The reference isn't that big...
        child = job.addChildJobFn(make_chromosome_graph_job, options, plan,
            chromosome_name, cores="1", memory="10G", disk="50G")
            
        contig_graph_ids.append(child.rv())
        
    return contig_graph_ids
            
        
            
            
def make_chromosome_graph_job(job, options, plan, chromosome_name):
    """
    Download the VCF and FASTAs for the given chromosome from the Toil job
    store, and use vg to make a graph for the contig with the given name in VCF
    space. Resulting graph will have contigs named in FASTA accession.version
    namespace.
    
    """
    
    # Find the VCF and index
    vcf_id = plan.get_vcf_id(chromosome_name)
    vcf_index_id = plan.get_vcf_index_id(chromosome_name)
    
    # Download them
    vcf_filename = os.path.join(job.fileStore.getLocalTempDir(), "vcf.vcf.gz")
    job.fileStore.readGlobalFile(vcf_id, vcf_filename)
    vcf_index_filename = vcf_filename + ".tbi"
    job.fileStore.readGlobalFile(vcf_index_id, vcf_index_filename)
    
    # Find the FASTAs
    fasta_ids = list(plan.for_each_primary_fasta())
    
    # Download the FASTAs
    fasta_filenames = []
    
    # FASTA files need to be unzipped
    for i, fasta_id in enumerate(fasta_ids):
        # Decide where to put the uncompressed FASTA.
        # Make sure it's in a place where we can put the index next to it
        fasta_filename = os.path.join(job.fileStore.getLocalTempDir(),
            "fasta{}.fa".format(i))
        RealtimeLogger.info("Downloading FASTA {}...".format(i))
        with job.fileStore.readGlobalFileStream(fasta_id) as compressed_stream:
            # Decompress each FASTA, if necessary
            decompressed = TransparentUnzip(compressed_stream)
            shutil.copyfileobj(decompressed, open(fasta_filename, "w"))
        
        fasta_filenames.append(fasta_filename)

    # Compose some VG arguments. We want to build just this one contig.
    vg_args = ["construct", "--region", chromosome_name, "--region-is-chrom",
        "--vcf", vcf_filename]
    
    for vcf_contig, fasta_contig in plan.get_name_translation().iteritems():
        # Add translations for all contig names.
        # TODO: Maybe we can just get by with our one, if present?
        # TODO: how do we know "=" isn't used in the names???
        vg_args += ["--rename", "{}={}".format(vcf_contig, fasta_contig)]
    
    for fasta_filename in fasta_filenames:
        # Add all the FASTAs. TODO: is vg going to explode if multiple VGs try
        # and index the FASTAs at once? And can it handle gzipped FASTAs?
        vg_args += ["--reference", fasta_filename]
    
    # Where will we put the graph?    
    graph_filename = job.fileStore.getLocalTempFile()
    
    RealtimeLogger.info("Running vg with: {}".format(vg_args))
    
    # Build the graph
    with open(graph_filename, "w") as graph_file:
        subprocess.check_call(["vg"] + vg_args, stdout=graph_file)
    
    # Upload the graph
    return job.fileStore.writeGlobalFile(graph_filename)
        
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    logging.info("Running on Toil from {}".format(toil.__file__))
    
    # Start up Toil
    with Toil(options) as toil_instance:
        
        if toil_instance.options.restart:
            # We're re-running. Grab the root job return value from restart
            graph_ids = toil_instance.restart()
        else:
            # Run from the top
        
            # Build the plan on the head node
            plan = create_plan(options.assembly_url, options.vcfs_url)
        
            # Import all the files from the plan. Now the plan will hold Toil IDs
            # for data files, and actual info for metadata files.
            plan.bake(lambda url: toil_instance.importFile(url))
    
            # Make a root job
            root_job = Job.wrapJobFn(main_job, options, plan,
                cores=1, memory="1G", disk="1G")
            
            # Run the root job
            graph_ids = toil_instance.start(root_job)
        
        for i, graph_id in enumerate(graph_ids):
            # Export all the graphs as VG files in arbitrary order
            toil_instance.exportFile(graph_id, "file:./part{}.vg".format(i))
        
    print("Toil workflow complete")
    return 0
    
def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        
        
        
        
        
        
        
        
        
        

