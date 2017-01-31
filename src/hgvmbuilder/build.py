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
    parser.add_argument("--hal_genome", action="append", default=[],
        help="extract the given genome from HAL files")
        
    # Contig conversion
    parser.add_argument("--add_chr", action="store_true",
        help="add chr prefix when converting from vcf to graph path names")
    
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
    
    # Find all the sequences we want
    sequences = []
    
    for genome in options.hal_genome:
        # Find all the sequences in the genome we want
        sequences_string = options.drunner.call([["halStats", hal_file, "--sequences", genome]], check_output=True)
        # Get a list of them
        sequences += sequences_string.strip().split(",")
    
    # TODO: we assume that these sequence names are unique in the HAL even
    # without the genome
    
    # Build a command to keep only the given paths in a vg file, and also chop
    mod_args = ["vg", "mod", "-", "--remove-non-path", "--chop", "100"]
    for sequence_name in sequences:
        # Retain each path in the correct genome
        mod_args.append("--retain-path")
        mod_args.append(sequence_name)
        RealtimeLogger.info("Retaining sequence {} from genomes {}".format(sequence_name, str(options.hal_genome)))
        
    with job.fileStore.writeGlobalFileStream() as (vg_handle, vg_id):
        # Make a vg and retain only the parts involved in the requested genome.
        # Also sort it so IDs will start at 1 again.
        # Save directly into the file store.
        options.drunner.call([["hal2vg", hal_file, "--chop", "1000000", "--inMemory", "--onlySequenceNames"],
            mod_args,
            ["vg", "ids", "-s", "-"]], outfile=vg_handle)
            
        RealtimeLogger.info("Created VG graph {}".format(vg_id))
            
        return vg_id
        
def xg_index_job(job, options, plan, vg_id):
    """
    Index the given VG graph into an XG file. Returns the ID of the XG file.
    """
    
    # Download the VG file
    vg_name = job.fileStore.readGlobalFile(vg_id)
    
    # Make an XG name
    xg_name = os.path.join(job.fileStore.getLocalTempDir(), "index.xg")
    
    # Set up args for the VG call
    vg_args = ["vg", "index", "-x", xg_name, vg_name]
    
    # Run the call
    options.drunner.call([vg_args])
    
    # Upload the file
    return job.fileStore.writeGlobalFile(xg_name)
    
def merge_vgs_job(job, options, plan, vg_ids):
    """
    Merge the given VG graphs into one VG graph. Returns the ID of the merged VG
    graph.
    """
    
    RealtimeLogger.info("Merging {} graphs".format(len(vg_ids)))
    
    if len(vg_ids) == 0:
        # Can't make a graph from no graphs
        raise RuntimeError("Cannot merge empty list of VG graphs!")
    elif len(vg_ids) == 1:
        # No work to do
        return vg_ids[0]
    else:
        # Actually need to merge
        # Grab all the graphs
        vg_names = [job.fileStore.readGlobalFile(vg_id) for vg_id in vg_ids]
        # Make a new filename for the merged graph
        merged_name = os.path.join(job.fileStore.getLocalTempDir(), "merged.vg")
        
        # Set up a command to merge all the graphs and put them in a shared ID
        # space
        vg_args = ["vg", "ids", "-j"] + vg_names
        
        with job.fileStore.writeGlobalFileStream() as (vg_handle, vg_id):
            # Make a merged vg
            options.drunner.call([vg_args], outfile=vg_handle)
            
        return vg_id
        
def hals_and_vgs_to_vg_job(job, options, plan, hal_ids, vg_ids):
    """
    Given some input HAL graphs to be converted to VG format, and some input VG
    graphs, combine all the graphs into one VG file with a consistent ID space.
    Return the ID of that VG file.
    """
    
    # First convert HALs to VG
    
    # Get promises for converting all the HALs into VGs    
    converted_vgs = []
    
    # Remember all the child jobs that make them
    hal_children = []
    
    for hal_id in hal_ids:
        # Make a child to convert each HAL to VG
        child = job.addChildJobFn(hal2vg_job, options, plan, hal_id,
            cores=1, memory="100G", disk="5G")
        converted_vgs.append(child.rv())
        hal_children.append(child)
        
    RealtimeLogger.info("Converting {} HALs".format(len(hal_children)))

    # Then make another child to merge all the graphs
    merge_child = job.addChildJobFn(merge_vgs_job, options, plan,
        converted_vgs + vg_ids, cores=1, memory="100G", disk="20G")
    for child in hal_children:
        # Make sure the merger comes after all the HAL imports
        child.addFollowOn(merge_child)
        
    return merge_child.rv()
    
def add_variants_job(job, options, plan, vg_id, vcf_ids):
    """
    Adds the VCF files with the given IDs into the VG file with the given ID
    using vg add. Returns the ID of a new vg file.
    
    """
    
    if len(vcf_ids) == 0:
        # Just spit back the graph unchanged
        return vg_id
    
    # Download all the VCFs
    vcf_filenames = []
    for i, vcf_id in enumerate(vcf_ids):
        # Pick a filename for each VCF
        vcf_filename = os.path.join(job.fileStore.getLocalTempDir(),
            "vcf{}.vcf.gz".format(i))
        
        # Download it
        job.fileStore.readGlobalFile(vcf_id, vcf_filename)
        # Download its index next to it
        job.fileStore.readGlobalFile(plan.get_index_id(vcf_id),
            vcf_filename + ".tbi")
            
        # Save the filename
        vcf_filenames.append(vcf_filename)
        
        RealtimeLogger.info("Downloaded VCF {} as {}".format(vcf_id,
            vcf_filename))
    
    # Download the input graph
    vg_filename = job.fileStore.readGlobalFile(vg_id)
    
    # Set up a command to add the VCFs to the graph
    vg_args = ["vg", "add", vg_filename]
    
    for vcf_filename in vcf_filenames:
        vg_args.append("-v")
        vg_args.append(vcf_filename)
        
        
    if not hasattr(options, 'add_chr') or options.add_chr:
        # Make sure to convert vcf names like "22" to UCSC-style "chr22".
        for base_name in plan.for_each_chromosome():
            vg_args.append("-n")
            vg_args.append("{}=chr{}".format(base_name, base_name))
    
    RealtimeLogger.info("Adding VCFs to vg graph...")
    
    with job.fileStore.writeGlobalFileStream() as (vg_handle, new_vg_id):
        # Stream new graph to the filestore
        options.drunner.call([vg_args], outfile=vg_handle)
    
    return new_vg_id
    
    
def main_job(job, options, plan):
    """
    Root Toil job. Returns a pair of lists of a VG file IDs and XG index file
    IDs.
    """
    
    # Make a child to convert the HALs and merge the VGs
    vg_job = job.addChildJobFn(hals_and_vgs_to_vg_job, options, plan,
        list(plan.for_each_hal()), list(plan.for_each_base_vg()),
        cores=1, memory="2G", disk="1G")
    
    # Add a followon to add in all the VCFs
    add_job = vg_job.addFollowOnJobFn(add_variants_job, options, plan,
        vg_job.rv(), [v for k, v in plan.for_each_vcf_id_by_chromosome()],
        cores=1, memory="100G", disk="100G")
    
    # Add a followon to index it
    xg_job = add_job.addFollowOnJobFn(xg_index_job, options, plan, add_job.rv(),
        cores=1, memory="100G", disk="20G")
    
    # Return the graph and the index
    return [add_job.rv()], [xg_job.rv()]
    
    
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
            graph_ids, index_ids = toil_instance.restart()
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
            graph_ids, index_ids = toil_instance.start(root_job)
        
        for i, graph_id in enumerate(graph_ids):
            # Export all the graphs as VG files in arbitrary order
            toil_instance.exportFile(graph_id, "file:./part{}.vg".format(i))
        for i, index_id in enumerate(index_ids):
            # Export all the graph indexes as XG files in the same order
            toil_instance.exportFile(index_id, "file:./part{}.xg".format(i))
        
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
        
        
        
        
        
        
        
        
        
        

