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
from .commandrunner import CommandRunner

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
    parser.add_argument("--chromosome", default=None,
        help="If specified, restrict to a single chromosome (like \"22\")")
    
    # Data sources
    parser.add_argument("--hal_url", default=None,
        help="start with the given HAL file and convert it to vg")
    parser.add_argument("--base_vg_url", default=None,
        help="start with this VG file instead of converting from HAL")
    parser.add_argument("--assembly_url", default=None,
        help="root of input assembly structure in GRC format")
    parser.add_argument("--vcfs_url", default=None,
        help="directory of VCFs per chromosome") 
        
    # HAL interpretation
    parser.add_argument("--hal_genome", default="Human",
        help="extract the given genome from HAL files")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def create_plan(assembly_url, vcfs_url, hal_url, base_vg_url):
    """
    Given an FTP or file url to the root of a GRC-format assembly_structure
    directory tree, and an FTP or file URL to a directory of chrXXX VCFs,
    produce a ReferencePlan describing that assembly and those VCFs.
    """
    
    # Make the plan
    plan = ReferencePlan()

    if assembly_url is not None:
        # Parse the assembly and populate the plan    
        grcparser.parse(plan, assembly_url)
    
    if vcfs_url is not None:
        # Parse the VCF directory and add the VCFs
        thousandgenomesparser.parse(plan, vcfs_url)
    
    if hal_url is not None:
        # We have a HAL file to convert to VG
        plan.add_hal(hal_url)
        
    if base_vg_url is not None:
        # We have a VG file already converted from HAL
        plan.add_base_vg(base_vg_url)
    
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
    
def concat_lists_job(job, lists):
    """
    Given a list of lists, concatenate them and return the result.
    
    Used to get a promise for the concatenation of promised lists.
    
    """
    
    RealtimeLogger.info("Concatenating {} lists".format(len(lists)))
    
    concatenated = []
    
    for l in lists:
        for item in l:
            # Put each item in each list into the big list
            concatenated.append(item)
            
    return concatenated
    
def hal2vg_job(job, options, plan, hal_id):
    """
    Convert a HAL file to a VG file and return the ID of the resulting vg graph.
    
    """
    
    # Download the HAL
    hal_file = job.fileStore.readGlobalFile(hal_id)
    
    RealtimeLogger.info("Converting HAL {}".format(hal_id))
    
    # Find all the sequences in the genome we want
    sequences_string = options.drunner.call([["halStats", hal_file, "--sequences", options.hal_genome]], check_output=True)
    # Get a list of them
    sequences = sequences_string.strip().split(",")
    
    # TODO: we assume that these sequence names are unique in the HAL even
    # without the genome
    
    # Build a command to keep only the given paths in a vg file, and also chop
    mod_args = ["vg", "mod", "-", "--remove-non-path", "--chop", "100"]
    for sequence_name in sequences:
        # Retain each path in the correct genome
        mod_args.append("--retain-path")
        mod_args.append(sequence_name)
        RealtimeLogger.info("Retaining sequence {} from genome {}".format(sequence_name, options.hal_genome))
        
    with job.fileStore.writeGlobalFileStream() as (vg_handle, vg_id):
        # Make a vg and retain only the parts involved in the requested genome.
        # Save directly into the file store.
        options.drunner.call([["hal2vg", hal_file, "--chop", "1000000", "--inMemory", "--onlySequenceNames"],
            mod_args], outfile=vg_handle)
            
        RealtimeLogger.info("Created VG graph {}".format(vg_id))
            
        return vg_id
    
    
def main_job(job, options, plan):
    """
    Root Toil job. Returns a list of file IDs for parts of the constructed
    graph.
    """
    
    # This list will hold lists of vg graphs, or promises of lists of vg graphs,
    # before any variants have been added.
    base_vgs = []
    
    # If we have base VGs, put them in the list
    base_vgs.append(list(plan.for_each_base_vg()))
    
    RealtimeLogger.info("Loading {} VGs".format(len(base_vgs[0])))

    # Get promises for converting all the HALs into VGs    
    hal_promises = []
    
    # Remember all the child jobs that make them
    hal_children = []
    
    for hal_id in plan.for_each_hal():
        # Make a child to convert each HAL to VG
        child = job.addChildJobFn(hal2vg_job, options, plan, hal_id,
            cores=1, memory="100G", disk="5G")
        hal_promises.append(child.rv())
        hal_children.append(child)
        
    RealtimeLogger.info("Converting {} HALs".format(len(hal_children)))

    # Add those to the list of vg graphs we have to start with
    base_vgs.append(hal_promises)

    # TODO: process the VG graphs and add in variants to them somehow.
    
    # Make a master list of graphs to return
    concat_job = job.addChildJobFn(concat_lists_job, base_vgs,
        cores=1, memory="1G", disk="1G")
    for before_child in hal_children:
        # Make sure HAL children run before this job that uses their RVs
        before_child.addFollowOn(concat_job)
    
    # Then return the list of all the VGs we have converted from HAL or taken in
    # to start with.
    return concat_job.rv()
        
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Cram a CommandRunner into the options, so our internal API looks like
    # toil-vg's.
    options.drunner = CommandRunner()
    
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Make sure VG is available
    options.drunner.call([["vg", "version"]])
    
    # And samtools
    # TODO: use vg to index FASTAs once, somehow.
    options.drunner.call([["samtools", "--version"]])
    
    # And hal2vg (which can't be made to succeed with no real arguments...)
    options.drunner.call([["which", "hal2vg"]])
    
    logging.info("Running on Toil from {}".format(toil.__file__))
    
    # Start up Toil
    with Toil(options) as toil_instance:
        
        if toil_instance.options.restart:
            # We're re-running. Grab the root job return value from restart
            graph_ids = toil_instance.restart()
        else:
            # Run from the top
        
            # Build the plan on the head node
            plan = create_plan(options.assembly_url, options.vcfs_url,
                options.hal_url, options.base_vg_url)
        
            # Import all the files from the plan. Now the plan will hold Toil
            # IDs for data files, and actual info for metadata files.
            plan.bake(lambda url: toil_instance.importFile(url))
        
            # Make a root job
            root_job = Job.wrapJobFn(main_job, options, plan,
                cores=1, memory="1G", disk="1G")
            
            # Run the root job and get one or more IDs for vg graphs, in a list
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
        
        
        
        
        
        
        
        
        
        

