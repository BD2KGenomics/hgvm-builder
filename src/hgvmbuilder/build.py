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
import re
import shutil
import subprocess
import sys

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
    
def prepare_fasta(job, fasta_id):
    """
    Given a job and a FASTA ID in its file store, download, decompress, scan,
    index, and re-upload the FASTA.
    
    Returns a tuple of the uncompressed FASTA file ID, the index file ID, and
    the list of accessions in the FASTA.
    
    Could be a Toil job or a function.
    """
    
    # Remember the contig names
    accessions = []
    
    RealtimeLogger.info("Preparing FASTA {}".format(fasta_id))
    
    # For each FASTA, drop it on disk somewhere
    filename = os.path.join(job.fileStore.getLocalTempDir(), "fasta.fa")
    with open(filename, "w") as out_handle:
        with job.fileStore.readGlobalFileStream(fasta_id) as in_handle:
            # Read from the file store, uncompress if needed, and write to disk.
            shutil.copyfileobj(TransparentUnzip(in_handle), out_handle)
         
    # Now index it
    subprocess.check_call(["samtools", "faidx", filename])   
    
    with open(filename + ".fai") as index_file:
        for line in index_file:
            # Crib the contig presence info from the index.
            accession = line.split("\t")[0]
            accessions.append(accession)
            RealtimeLogger.info("Found contig {}".format(accession))
    
    # Upload and return
    return (job.fileStore.writeGlobalFile(filename),
        job.fileStore.writeGlobalFile(filename + ".fai"), accessions)
    
    
    
def main_job(job, options, plan):
    """
    Root Toil job. Execute the plan. Returns a list of graph file IDs to export.
    """
    
    # This holds uncompressed FASTA ID, index ID, and contig list tuples from
    # prepared primary FASTAs
    prepared_primary = []
    # And this holds the same things for alt FASTAs
    prepared_alt = []
    
    # This holds child jobs that process FASTAs
    fasta_children = []
    
    for fasta_id in plan.for_each_primary_fasta_id():
        # Prepare each primary FASTA
        child = job.addChildJobFn(prepare_fasta, fasta_id,
            cores="1", memory="4G", disk="5G")
        fasta_children.append(child)
        prepared_primary.append(child.rv())
        
    
    for fasta_id in plan.for_each_alt_fasta_id():
        # Prepare each alt FASTA
        child = job.addChildJobFn(prepare_fasta, fasta_id,
            cores="1", memory="4G", disk="5G")
        fasta_children.append(child)
        prepared_alt.append(child.rv())
            
    # Add a child that will use those FASTAs, and return what it returns
    last_child = job.addChildJobFn(find_fastas_and_run_builds_job, options,
        plan, prepared_primary, prepared_alt, cores="1", memory="4G", disk="1G")
        
    for child in fasta_children:
        # Make sure all the other children are done before this last one can
        # start
        child.addFollowOn(last_child)
    
    # Return what that last child returns
    return last_child.rv()
    


def find_fastas_and_run_builds_job(job, options, plan, prepared_primary, prepared_alt):
    """
    Given the list of results from preparing all the primary FASTAs, and the
    list of results from preparing all the alt FASTAs, run the graph build.
    Returns a list of file IDs for vg graphs to export.
    """
    
    # This maps from primary accession to uncompressed fasta ID and index ID.
    primary_to_fasta = {}
    
    # This maps from alt accession to uncompressed fasta ID and index ID.
    alt_to_fasta = {}
    
    for uncompressed_id, index_id, contigs in prepared_primary:
        # Look at all the prepared primary FASTAs
        for contig in contigs:
            # Every contig found in that FASTA points to the uncompressed FASTA
            # and its index.
            primary_to_fasta[contig] = (uncompressed_id, index_id)
            
    
    for uncompressed_id, index_id, contigs in prepared_alt:
        # Look at all the prepared alt FASTAs
        for contig in contigs:
            # Every contig found in that FASTA points to the uncompressed FASTA
            # and its index.
            alt_to_fasta[contig] = (uncompressed_id, index_id)
           
    # Grab the dict of VCF IDs by chromosome names
    # TODO: stop converting back and forth
    vcf_ids_by_chromosome = dict(plan.for_each_vcf_id_by_chromosome())
            
    # Build a list of VG graph promises
    contig_graph_ids = []
    
    RealtimeLogger.info("Building graph from {} primary contigs and {} "
        "VCFs".format(len(primary_to_fasta), len(vcf_ids_by_chromosome)))
            
    for (accession, (fasta_id, fasta_index_id)) in primary_to_fasta.iteritems():
        # For each primary accession, grab its FASTA ID and index ID
        
        # Get the chromosome name (or just the accession again if it doesn't
        # have one)
        chromosome_name = plan.accession_to_chromosome_name(accession)
        
        # Grab its VCF ID and VCF index ID, if any
        vcf_id = None
        vcf_index_id = None
        
        if vcf_ids_by_chromosome.has_key(chromosome_name):
            # We catually have these files
            vcf_id = vcf_ids_by_chromosome[chromosome_name]
            vcf_index_id = plan.get_index_id(vcf_id)
        
        # Dispatch a child to build it with those four files and the contig name
        # (or accession if there is no contig name).
        child = job.addChildJobFn(make_chromosome_graph_job, options, plan,
            chromosome_name, fasta_id, fasta_index_id, vcf_id, vcf_index_id,
            cores="1", memory="10G", disk="50G")
            
        contig_graph_ids.append(child.rv())
        
    return contig_graph_ids
            
        
            
            
def make_chromosome_graph_job(job, options, plan, chromosome_name, fasta_id,
    fasta_index_id, vcf_id, vcf_index_id):
    """
    Download the FASTA and its index, and the VCF and its index if given, and
    use vg to make a graph for the contig with the given name in VCF space.
    Resulting graph will have contigs named in FASTA accession.version
    namespace.
    
    """
    
    # Compose some VG arguments. We want to build just this one contig.
    vg_args = ["construct", "--region", chromosome_name, "--region-is-chrom"]
    
    if vcf_id is not None:
        # We have a VCF and an index. Download them.
        vcf_filename = os.path.join(job.fileStore.getLocalTempDir(), "vcf.vcf.gz")
        job.fileStore.readGlobalFile(vcf_id, vcf_filename)
        vcf_index_filename = vcf_filename + ".tbi"
        job.fileStore.readGlobalFile(vcf_index_id, vcf_index_filename)
        
        vg_args += ["--vcf", vcf_filename]
        
    # Download the FASTA and FASTA index
    fasta_filename = os.path.join(job.fileStore.getLocalTempDir(), "fasta.fa")
    job.fileStore.readGlobalFile(fasta_id, fasta_filename)
    fasta_index_filename = fasta_filename + ".fai"
    job.fileStore.readGlobalFile(fasta_index_id, fasta_index_filename)
    vg_args += ["--reference", fasta_filename]

    # Add the necessary rename for this chromosome to its accession
    vg_args += ["--rename", "{}={}".format(chromosome_name,
        plan.chromosome_name_to_accession(chromosome_name))]

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
    
    # Make sure VG is available
    subprocess.check_call(["vg", "version"])
    
    # And samtools
    # TODO: use vg to index FASTAs once, somehow.
    subprocess.check_call(["samtools", "--version"])
    
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
        
        
        
        
        
        
        
        
        
        

