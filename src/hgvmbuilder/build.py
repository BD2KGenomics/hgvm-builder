#!/usr/bin/env python
# hgvm-builder build.py: Command-line tool to build a Human Genome Variation Map
"""

Takes a GRC-format assembly strucutre, and a VCF, and produces a vg graph of the
assembly (including novel and fix patches, alternate loci, and
unlocalized/poorly localized scaffonds) including the variation from the VCF.

Parallelizes using the Toil system.

Requires the "vg" binary to be available on the PATH on all nodes.

"""

import argparse
import logging
import os
import os.path
import re
import shutil
import subprocess
import sys
import collections

import toil
from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from .plan import ReferencePlan
from .evaluation import EvaluationPlan
from . import grcparser
from . import thousandgenomesparser
from .transparentunzip import TransparentUnzip
from .commandrunner import CommandRunner
from .directory import Directory

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
    # For construction
    parser.add_argument("--hal_url", default=None,
        help="start with the given HAL file and convert it to vg")
    parser.add_argument("--base_vg_url", default=None,
        help="start with this VG file instead of converting from HAL")
    parser.add_argument("--assembly_url", default=None,
        help="root of input assembly structure in GRC format")
    parser.add_argument("--vcfs_url", action="append", default=[],
        help="directory of VCFs per chromosome (may repeat)")
    # For evaluation
    parser.add_argument("--control_graph_url", default=None,
        help="use the given vg file as a control")
    parser.add_argument("--control_graph_xg_url", default=None,
        help="use the given xg index for the control graph")
    parser.add_argument("--control_graph_gcsa_url", default=None,
        help="use the given gcsa/lcp index for the control graph")
    parser.add_argument("--sample_fastq_url", action="append", default=[],
        help="FASTQ for one end of evaluation reads")
    parser.add_argument("--eval_sequences_url", default=None,
        help="sequences to align to the sample graph, one sequence per line")
        
        
    # HAL interpretation
    parser.add_argument("--hal_genome", action="append", default=[],
        help="extract the given genome from HAL files")
        
    # Contig conversion
    parser.add_argument("--add_chr", action="store_true",
        help="add chr prefix when converting from vcf to graph path names")
        
    # Output
    parser.add_argument("out_dir",
        help="directory to place the constructed graph parts in")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def create_plan(assembly_url, vcfs_urls, hal_url, base_vg_url):
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
    
    for vcfs_url in vcfs_urls:
        # Parse each VCF directory and add the VCFs
        thousandgenomesparser.parse(plan, vcfs_url)
    
    if hal_url is not None:
        # We have a HAL file to convert to VG
        plan.add_hal(hal_url)
        
    if base_vg_url is not None:
        # We have a VG file already converted from HAL
        plan.add_base_vg(base_vg_url)
    
    # Return the completed plan
    return plan
    
def create_eval_plan(control_url, xg_url, gcsa_url, sample_fq_urls,
    sequences_url):
    """
    Create an EvaluationPlan to compare againbst the given control graph (with
    the given xg, and GCSA/LCP indexes), by variant calling with the given
    FASTQs and then realigning the given sequences to the sample graph.
    """
    
    eval_plan = EvaluationPlan()
    
    if control_url is not None:
        # Grab the control graph
        eval_plan.set_control_graph(control_url)
    
    if xg_url is not None:
        # And its XG index
        eval_plan.set_control_graph_xg(xg_url)
   
    if gcsa_url is not None:
        # And its GCSA/LCP index. Assume the lcp is named right
        eval_plan.set_control_graph_gcsa(gcsa_url, gcsa_url + ".lcp")
        
    for url in sample_fq_urls:
        # And all the FASTQs
        eval_plan.add_fastq(url)
        
    if sequences_url is not None:
        # And the sequences that will be used to evaluate the called graphs
        eval_plan.set_eval_sequences(sequences_url)
        
    return eval_plan
   
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
    
    # We need to want a HAL genome to use a HAL
    assert(len(options.hal_genome) > 0)
    
    # Download the HAL
    hal_file = job.fileStore.readGlobalFile(hal_id)
    
    RealtimeLogger.info("Converting HAL {}".format(hal_id))
    
    # Find all the sequences we want
    sequences = []
    
    # Convert into a local file
    vg_filename = os.path.join(job.fileStore.getLocalTempDir(), "from_hal.vg")
    
    # Pull out only the requested genomes, and sort so IDs start at 1 and are
    # nice.
    options.drunner.call([["hal2vg", hal_file,
        "--inMemory",
        "--onlySequenceNames",
        "--targetGenomes", ",".join(options.hal_genome),
        "--refGenome", options.hal_genome[0],
        "--noAncestors"],
        ["vg", "mod", "-", "--chop", "100"],
        ["vg", "ids", "-s", "-"]], outfile=open(vg_filename, "w"))
        
    # Validate the resulting graph
    options.drunner.call([["vg", "validate", vg_filename]])
    
    # Upload it
    vg_id = job.fileStore.writeGlobalFile(vg_filename)
        
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
    
def gcsa_index_job(job, options, plan, vg_id):
    """
    Index the given VG graph into a GCSA2/LCP file pair. Returns the IDs of the
    .gcsa and .lcp files in a pair.
    """
    
    # Download the VG file
    vg_name = job.fileStore.readGlobalFile(vg_id)
    
    # Make a GCSA name
    gcsa_name = os.path.join(job.fileStore.getLocalTempDir(), "index.gcsa")
    # Find the LCP name
    lcp_name = gcsa_name + ".lcp"
    
    # Set up args for the VG call
    # TODO: don't hardcode options
    vg_args = ["vg", "index", "-g", gcsa_name,
        "-k", "8", "-e", "2", "-X", "2", "--size-limit", "500", vg_name]
    
    # Run the call
    options.drunner.call([vg_args])
    
    # Upload the files
    return (job.fileStore.writeGlobalFile(gcsa_name),
        job.fileStore.writeGlobalFile(lcp_name))
    
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
        # Actually need to merge.
        
        # Grab all the graphs. Don't cache because we'll need to modify their
        # IDs in place.
        vg_names = [job.fileStore.readGlobalFile(vg_id, cache=False)
            for vg_id in vg_ids]
        # Make a new filename for the merged graph
        merged_name = os.path.join(job.fileStore.getLocalTempDir(), "merged.vg")
        
        for vg_filename in vg_names:
            # Make sure all the vg files are writeable
            # Cache=False might be supposed to do it but we have to make sure.
            os.chmod(vg_filename, 0644)
        
        # Set up a command to merge all the graphs and put them in a shared ID
        # space
        vg_args = ["vg", "ids", "-j"] + vg_names
        
        # It modifies the files in place
        options.drunner.call([vg_args])
        
        cat_args = ["cat"] + vg_names
        
        with job.fileStore.writeGlobalFileStream() as (vg_handle, vg_id):
            # Make a merged vg with cat
            options.drunner.call([cat_args], outfile=vg_handle)
            
        return vg_id
        
def explode_vg_job(job, options, plan, vg_id):
    """
    Explode the given vg file into connected component subgraphs.
    
    Returns a dict from subgraph ID to a list of path names present in the
    subgraph.
    """
    
    RealtimeLogger.info("Exploding graph {}".format(vg_id))
    
    # Download the VG file
    vg_name = job.fileStore.readGlobalFile(vg_id)
    
    # Set a directory for parts
    part_dir = job.fileStore.getLocalTempDir() + "/parts"
    
    # Reserve a report file
    report_filename = job.fileStore.getLocalTempDir() + "/report.tsv"
    
    # Set up args for the VG explode call
    vg_args = ["vg", "explode", vg_name, part_dir]
    
    # Run the call and capture the output
    report = options.drunner.call([vg_args], outfile=open(report_filename, "w"))
    
    RealtimeLogger.info("Exploded and put report in {}".format(report_filename))
    
    # We'll populate this with the paths that live in each part
    paths_by_part = {}
    
    for line in open(report_filename):
        if line == "":
            # Skip blank lines
            continue
        # Split each line on tabs
        parts = line.split("\t")
        
        # Validate the part
        options.drunner.call([["vg", "validate", parts[0]]])
        
        # Save the file (first entry) to the filestore
        part_id = job.fileStore.writeGlobalFile(parts[0])
        
        # Remember that all of its paths belong to it
        paths_by_part[part_id] = parts[1:]
        
        RealtimeLogger.info("Placed {} paths in part {}".format(len(parts) - 1, parts[0]))
        
    assert len(paths_by_part) > 0
        
    return paths_by_part
        
def add_variants_to_parts_job(job, options, plan, paths_by_part):
    """
    Given a dict from vg file ID to the paths in that VG file, work out the VCF
    files for each VG file and run jobs to add their variants. Return a list of
    vg IDs for the graphs with variants added.
    
    """
    
    # This will hold VCF IDs by the path name they apply to
    vcfs_by_path = collections.defaultdict(list)
    
    for chromosome, vcf_id in plan.for_each_vcf_id_by_chromosome():
        # Build a dict from path name to relevant VCFs
        
        if not hasattr(options, 'add_chr') or options.add_chr:
            # Make sure to convert vcf names like "22" to UCSC-style "chr22".
            chromosome = "chr" + chromosome
            
        vcfs_by_path[chromosome].append(vcf_id)
    
    # Make a list of promises for updated vg graph file IDs
    to_return = []
    
    for vg_id, paths in paths_by_part.iteritems():
        # For each VG graph
        
        # Grab the applicable VCFs, flattening across all paths in the VG
        vcfs = list({vcf_id for path in paths for vcf_id in vcfs_by_path[path]})
        
        # Add those VCFs to this VG
        child = job.addChildJobFn(add_variants_job, options, plan, vg_id,
            vcfs, cores=2, memory="100G", disk="50G")
            
        to_return.append(child.rv())
        
    return to_return
        
        
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
    
    if len(options.hal_genome) > 0:
        # We want something from the HALs
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
        # We need to download each VCF
        
        if plan.get_index_id(vcf_id) is not None:
            # This VCF is indexed, so it must be gzipped
            
            # Download the VCF under a .gz name
            vcf_filename = os.path.join(job.fileStore.getLocalTempDir(),
            "vcf{}.vcf.gz".format(i))
            job.fileStore.readGlobalFile(vcf_id, vcf_filename)
            
            # Download its index next to it
            job.fileStore.readGlobalFile(plan.get_index_id(vcf_id),
                vcf_filename + ".tbi")
            
            # Make sure the index is newer than the VCF
            os.utime(vcf_filename + ".tbi", None)
            
        else:
            # If it's not indexed, assume it's not compressed either.
            vcf_filename = os.path.join(job.fileStore.getLocalTempDir(),
            "vcf{}.vcf".format(i))
            job.fileStore.readGlobalFile(vcf_id, vcf_filename)
            
        # Save the filename
        vcf_filenames.append(vcf_filename)
        
        RealtimeLogger.info("Downloaded VCF {} as {}".format(vcf_id,
            vcf_filename))
    
    # Download the input graph
    vg_filename = job.fileStore.readGlobalFile(vg_id)
    
    # Validate it
    options.drunner.call([["vg", "validate", vg_filename]])
    
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
    
    RealtimeLogger.info("Adding {} VCFs to vg graph...".format(
        len(vcf_filenames)))
    
    with job.fileStore.writeGlobalFileStream() as (vg_handle, new_vg_id):
        # Stream new graph to the filestore
        options.drunner.call([vg_args], outfile=vg_handle)
    
    return new_vg_id
    
def hgvm_package_job(job, options, plan, vg_id, xg_id, gcsa_pair):
    """
    Given a VG file ID, and XG file ID, and a pair of GCSA and LCP file IDs,
    produce a Directory with "hgvm.vg", "hgvm.xg", "hgvm.gcsa", and
    "hgvm.gcsa.lcp" filled in.
    
    Runs as a Toil job so it will be able to unpack the pair.
    
    """
    
    # Package up the file IDs into a Directory
    hgvm = Directory()
    hgvm.add("hgvm.vg", vg_id)
    hgvm.add("hgvm.xg", xg_id)
    hgvm.add("hgvm.gcsa", gcsa_pair[0])
    hgvm.add("hgvm.gcsa.lcp", gcsa_pair[1])
    
    return hgvm
    
    
    
def hgvm_build_job(job, options, plan):
    """
    Toil job for building the HGVM. Returns a Directory object holding
    "hgvm.vg", "hgvm.xg", "hgvm.gcsa", and "hgvm.gcsa.lcp".
    
    """
    
    RealtimeLogger.info("Building HGVM")
    
    # Make a child to convert the HALs and merge the VGs
    vg_job = job.addChildJobFn(hals_and_vgs_to_vg_job, options, plan,
        list(plan.for_each_hal()), list(plan.for_each_base_vg()),
        cores=1, memory="2G", disk="1G")
        
    # Then split the VG up into connected components again
    explode_job = vg_job.addFollowOnJobFn(explode_vg_job, options, plan,
        vg_job.rv(),
        cores=1, memory="100G", disk="100G")
        
    # And add all the appropriate VCFs to each
    add_job = explode_job.addFollowOnJobFn(add_variants_to_parts_job, options,
        plan, explode_job.rv(),
        cores=1, memory="4G", disk="4G")
        
    # And then merge those VCFs together
    merge_job = add_job.addFollowOnJobFn(merge_vgs_job, options, plan,
        add_job.rv(),
        cores=1, memory="100G", disk="100G")
    
    # Add a followon to XG-index it
    xg_job = merge_job.addFollowOnJobFn(xg_index_job, options, plan,
        merge_job.rv(),
        cores=1, memory="100G", disk="20G")
    
    # And another to GCSA-index it
    gcsa_job = merge_job.addFollowOnJobFn(gcsa_index_job, options, plan,
        merge_job.rv(),
        cores=16, memory="100G", disk="500G")
    
    # And a job to package it all up, depending on the XG and GCSA.
    directory_job = xg_job.addFollowOnJobFn(hgvm_package_job, options, plan,
        merge_job.rv(), xg_job.rv(), gcsa_job.rv(),
        cores=1, memory="2G", disk="1G")
    gcsa_job.addFollowOn(directory_job)
    
    # Return the Directory with the packaged HGVM
    return directory_job.rv()
    
def hgvm_eval_job(job, options, eval_plan, hgvm):
    """
    Toil job for evaluating the HGVM. Takes a packaged HGVM Directory, and a
    plan for evaluating it.
    
    Returns the evaluation results as a Directory.
    
    """
    
    RealtimeLogger.info("Evaluating HGVM")
    
    results = Directory()
    
    # TODO: implement
    
    return results
    
def main_job(job, options, plan, eval_plan):
    """
    Root job of the Toil workflow. Build and then evaluate an HGVM.
    
    Returns a pair of the packaged HGVM Directory and the evaluation results
    Directory.
    
    """
    
    RealtimeLogger.info("Main job starting")
    
    # Build the HGVM
    build_job = job.addChildJobFn(hgvm_build_job, options, plan,
        cores=1, memory="2G", disk="1G")
    
    # Then evaluate it
    eval_job = build_job.addFollowOnJobFn(hgvm_eval_job, options, eval_plan,
        build_job.rv(),
        cores=1, memory="2G", disk="1G")
        
    # Return the pair of created Directories
    return (build_job.rv(), eval_job.rv())
    
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
                
            # Also build the evaluation plan on the head node
            eval_plan = create_eval_plan(options.control_graph_url,
                options.control_graph_xg_url, options.control_graph_gcsa_url,
                options.sample_fastq_url, options.eval_sequences_url)
        
            # Import all the files from the plans. Now the plans will hold Toil
            # IDs for data files, and actual info for metadata files.
            plan.bake(lambda url: toil_instance.importFile(url))
            eval_plan.bake(lambda url: toil_instance.importFile(url))
        
            # Make a root job
            root_job = Job.wrapJobFn(main_job, options, plan, eval_plan,
                cores=1, memory="1G", disk="1G")
            
            # Run the root job and get the packaged HGVM and its evaluation as
            # Directories.
            hgvm_directory, eval_directory = toil_instance.start(root_job)
            
        try:
            # Make sure the output directory exists
            os.makedirs(options.out_dir)
        except OSError:
            # It may exist already
            pass
            
        # Export the HGVM
        hgvm_directory.export(toil_instance, "file:{}".format(options.out_dir))
        # And stick its evaluations inside it for now
        eval_directory.export(toil_instance,
            "file:{}/eval".format(options.out_dir))
        
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
        
        
        
        
        
        
        
        
        
        

