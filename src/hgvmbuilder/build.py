#!/usr/bin/env python
# hgvm-builder build.py: Command-line tool to build a Human Genome Variation Map
"""

Takes a GRC-format assembly strucutre, and a VCF, and produces a vg graph of the
assembly (including novel and fix patches, alternate loci, and
unlocalized/poorly localized scaffonds) including the variation from the VCF.

Parallelizes using the Toil system.

Requires the "vg", "samtools", and "hal2vg" binaries to be available on the PATH
on all nodes, or to use Docker or Singularity.

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
import itertools
import more_itertools
import shutil
import datetime
import uuid
import random

import tsv
import intervaltree
import vcf

import toil
import toil.version
from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from .plan import ReferencePlan
from .evaluation import EvaluationPlan, SVRecallPlan
from . import grcparser
from . import thousandgenomesparser
from .transparentunzip import TransparentUnzip
from .commandrunner import CommandRunner
from .directory import Directory
from .toilpromise import ToilPromise
from .version import version
from . import toilvgfacade
from .manifest import Manifest

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
    
    # Add all the toil-vg options
    toilvgfacade.add_options(parser)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # Add a version
    parser.add_argument("--version", action="version",
        version="%(prog)s {}\nToil {}: {}\nPyVCF: {}".format(version,
            toil.version.version, toil.__file__, vcf.__file__),
        help="print version and exit")
    
    # General options
    parser.add_argument("--chromosome", default=None,
        help="If specified, restrict to a single chromosome (like \"22\")")
    
    # Data sources
    # For construction
    construction_group = parser.add_argument_group("Construction")
    construction_group.add_argument("--hal_url", default=None,
        help="start with the given HAL file and convert it to vg")
    construction_group.add_argument("--base_vg_url", default=None,
        help="start with this VG file instead of converting from HAL")
    construction_group.add_argument("--assembly_url", default=None,
        help="root of input assembly structure in GRC format")
    construction_group.add_argument("--vcfs_url", action="append", default=[],
        help="directory of VCFs per chromosome (may repeat)")
    construction_group.add_argument("--vcf_url", action="append", default=[],
        help="single VCF for a particular chromosome")
    construction_group.add_argument("--vcf_contig", action="append", default=[],
        help="contig/chromosome for the corresponding VCF")
    # HAL interpretation
    construction_group.add_argument("--hal_genome", action="append", default=[],
        help="extract the given genome from HAL files")
    # Contig conversion
    construction_group.add_argument("--add_chr", action="store_true",
        help="add chr prefix when converting from vcf to graph path names")
        
    # We can just skip construction proper and import an HGVM pre-packaged from
    # a directory.
    construction_group.add_argument("--hgvm_url", default=None,
        help="just import the given pre-built HGVM for evaluation")
        
    # For realignment evaluation
    realignment_group = parser.add_argument_group("Realignment evaluation")
    realignment_group.add_argument("--control_graph_hgvm_url", default=None,
        help="use the given hgvm directory as a control")
    realignment_group.add_argument("--control_graph_url", default=None,
        help="use the given vg file as a control")
    realignment_group.add_argument("--control_graph_xg_url", default=None,
        help="use the given xg index for the control graph")
    realignment_group.add_argument("--control_graph_gcsa_url", default=None,
        help="use the given gcsa/lcp index for the control graph")
    realignment_group.add_argument("--sample_fastq_url", action="append",
        default=[],
        help="FASTQ for one end of realignment evaluation reads")
    realignment_group.add_argument("--eval_sequences_url", default=None,
        help="sequences to align to the sample graph, one sequence per line")
        
    # For SV calling and comparison against truth evaluaton
    sv_group = parser.add_argument_group("Structural variant evaluation")
    sv_group.add_argument("--sv_sample_name", action="append", default=[],
        help="evaluate recall of this sample's SVs")
    
    # We need reads and we will get them form only one place
    sv_read_source = sv_group.add_mutually_exclusive_group()
    
    # This list is going to hold FASTQs for multiple samples; all samples must
    # be either paired-end or single-end so we can match them up.
    sv_read_source.add_argument("--sv_sample_fastq_url", action="append",
        default=[],
        help="FASTQ for one end of SV evaluation reads for one sample")
    # Instead of FASTQs we can use pre-aligned GAMs. But we can only use one or
    # the other.
    sv_read_source.add_argument("--sv_sample_gam_url", action="append",
        default=[],
        help="GAM for pre-aligned evaluation reads instead of FASTQs")
    # Instead of either of those, we could also use VCFs instead.
    sv_read_source.add_argument("--sv_sample_vcf_url", action="append",
        default=[],
        help="VCF for pre-computed sample genotypes for recall analysis")
        
    # Output
    parser.add_argument("out_url",
        help="file: or other Toil-supported URL to place results in")
    parser.add_argument("--dump_hgvm", default=None, type=os.path.abspath,
        help="dump the HGVM build to the given directory")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def create_plan(assembly_url, vcfs_urls, vcf_urls, vcf_contigs, hal_url,
    base_vg_url, hgvm_url):
    """
    Given an FTP or file url to the root of a GRC-format assembly_structure
    directory tree, a list of FTP or file URLs to directories of chrXXX VCFs, a
    list of URLs to individual VCFs to include, a list of contigs to assign
    those VCFs to, a HAL URL to import, and a pre-built HGVM directory URL to
    import, produce a ReferencePlan describing the HGVM we want.
    
    Generally not all of those will be in use at the same time.
    """
    
    # Make the plan
    plan = ReferencePlan()

    if assembly_url is not None:
        # Parse the assembly and populate the plan    
        grcparser.parse(plan, assembly_url)
    
    for vcfs_url in vcfs_urls:
        # Parse each VCF directory and add the VCFs
        thousandgenomesparser.parse(plan, vcfs_url)
        
    if len(vcf_contigs) != len(vcf_urls):
        # We need these to be in correspondence
        raise RuntimeError("Must use the same number of VCFs as contigs")
        
    for vcf_url, vcf_contig in itertools.izip(vcf_urls, vcf_contigs):
        # Add each single-contig VCF for its contig
        plan.add_variants(vcf_contig, vcf_url)
        
        if vcf_url.lower().endswith(".gz"):
            # Also add the index that we assume exists
            plan.add_variants_index(vcf_url + ".tbi")
    
    if hal_url is not None:
        # We have a HAL file to convert to VG
        plan.add_hal(hal_url)
        
    if base_vg_url is not None:
        # We have a VG file already converted from HAL
        plan.add_base_vg(base_vg_url)
        
    if hgvm_url is not None:
        # We just want to use this HGVM we have prepared earlier.
        plan.set_hgvm(hgvm_url)
    
    
    # Return the completed plan
    return plan
    
def create_eval_plan(control_hgvm_url, control_url, xg_url, gcsa_url,
    sample_fq_urls, sequences_url):
    """
    Create an EvaluationPlan to compare againbst the given control graph (with
    the given xg, and GCSA/LCP indexes), by variant calling with the given
    FASTQs and then realigning the given sequences to the sample graph.
    """
    
    eval_plan = EvaluationPlan()
    
    if control_hgvm_url is not None:
        # Grab the control hgvm
        eval_plan.set_control_hgvm(control_hgvm_url)
        
        # We can only use one of these sets. TODO: push this requirement into
        # argparse.
        assert(control_url is None)
        assert(xg_url is None)
        assert(gcsa_url is None)
    else:
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
    
def create_recall_plan(sample_fastq_urls, gam_urls, vcf_urls, sample_names):
    """
    Create an SVRecallPlan to evaluate structural variant recall for the given
    sample, using the given FASTQs or pre-aligned GAM.
    """
    
    recall_plan = SVRecallPlan()
    
    for sample_name in sample_names:
        # Add a sample by this name
        recall_plan.add_sample_name(sample_name)
        
    if len(sample_fastq_urls) == len(sample_names):
        # Unpaired FASTQ (or no FASTQs)
        for sample_name, sample_fastq_url in itertools.izip(sample_names,
            sample_fastq_urls):
            # Operate on one FASTQ per sample, and associate it with the sample
            recall_plan.add_sample_fastq(sample_name, sample_fastq_url)
    elif len(sample_fastq_urls) == 2 * len(sample_names):
        for sample_name, (fastq1, fastq2) in itertools.izip(sample_names,
            more_itertools.grouper(2, sample_fastq_urls)):
            # For every sample name and the corresponding pair of FASTQs
            # Add both the FASTQs for the sample
            recall_plan.add_sample_fastq(sample_name, fastq1)
            recall_plan.add_sample_fastq(sample_name, fastq2)
    elif len(sample_fastq_urls) != 0:
        # Not sure what the user thinks we're going to be able to do with this
        # weird mismatched set of FASTQs.
        raise RuntimeError("Cannot match up {} FASTQs with {} samples; "
            "provide either all single-end or all paired FASTQs.".format(
            len(sample_fastq_urls), len(sample_names)))
    
    if len(gam_urls) == len(sample_names):
        # We have one GAM per sample instead of FASTQs
        for sample_name, gam_url in itertools.izip(sample_names, gam_urls):
            # Operate on one GAM per sample, and associate it with the sample
            recall_plan.add_sample_gam(sample_name, gam_url)
    elif len(gam_urls) != 0:
        raise RuntimeError("Cannot match up {} GAMs with {} samples.".format(
            len(gam_urls), len(sample_names)))
            
    if len(vcf_urls) == len(sample_names):
        # We have one VCF per sample instead of FASTQs or GAMs
        for sample_name, vcf_url in itertools.izip(sample_names, vcf_urls):
            # Operate on one VCF per sample, and associate it with the sample
            recall_plan.add_sample_vcf(sample_name, vcf_url)
    elif len(vcf_urls) != 0:
        raise RuntimeError("Cannot match up {} VCFs with {} samples.".format(
            len(vcf_urls), len(sample_names)))
            
    if (len(gam_urls) == 0 and len(sample_fastq_urls) == 0 and
        len(vcf_urls) == 0 and len(sample_names) > 0):
        # We have no reads or calls for our samples!
        raise RuntimeError("Cannot do recall evaluation on any samples "
            "without reads or calls for them")
        
        
    return recall_plan
    
   
def prepare_fasta(job, options, fasta_id):
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
    
    # Make a work_dir we can mount in Docker
    work_dir = job.fileStore.getLocalTempDir()
    
    # For each FASTA, drop it on disk somewhere
    filename = os.path.join(work_dir, "fasta.fa")
    with open(filename, "w") as out_handle:
        with job.fileStore.readGlobalFileStream(fasta_id) as in_handle:
            # Read from the file store, uncompress if needed, and write to disk.
            shutil.copyfileobj(TransparentUnzip(in_handle), out_handle)
         
    # Now index it
    options.drunner.call([["samtools", "faidx", "fasta.fa"]], work_dir=work_dir)   
    
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
    
    # Define a work_dir for Docker purposes.
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the HAL
    hal_file = "halfile.hal"
    job.fileStore.readGlobalFile(hal_id, os.path.join(work_dir, hal_file))
    
    RealtimeLogger.info("Converting HAL {}".format(hal_id))
    
    # Find all the sequences we want
    sequences = []
    
    # Convert into a local file.
    vg_filename = "from_hal.vg"
    
    # Pull out only the requested genomes, and sort so IDs start at 1 and are
    # nice.
    options.drunner.call(job, [["hal2vg", hal_file,
        "--inMemory",
        "--onlySequenceNames",
        "--targetGenomes", ",".join(options.hal_genome),
        "--refGenome", options.hal_genome[0],
        "--noAncestors"],
        ["vg", "mod", "-", "--chop", "100"],
        ["vg", "ids", "-s", "-"]], outfile=open(
        os.path.join(work_dir, vg_filename), "w"), work_dir=work_dir)
        
    # Validate the resulting graph
    options.drunner.call(job, [["vg", "validate", vg_filename]],
        work_dir=work_dir)
    
    # Upload it
    vg_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, vg_filename))
        
    RealtimeLogger.info("Created VG graph {}".format(vg_id))
    return vg_id
        
def merge_vgs_job(job, options, vg_ids):
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
        
        work_dir = job.fileStore.getLocalTempDir()
        
        # Make up unique names for all the graphs we need
        vg_names = ["graph{}.vg".format(i) for i in xrange(len(vg_ids))]
        
        for name, file_id in itertools.izip(vg_names, vg_ids):
            # Grab all the graphs. Don't cache because we'll need to modify
            # their IDs in place.
            job.fileStore.readGlobalFile(file_id,
                userPath=os.path.join(work_dir, name), cache=False)

        # Make a new filename for the merged graph
        merged_name = "merged.vg"
        
        for vg_filename in vg_names:
            # Make sure all the vg files are writeable
            # Cache=False might be supposed to do it but we have to make sure.
            os.chmod(os.path.join(work_dir, vg_filename), 0644)
        
        # Set up a command to merge all the graphs and put them in a shared ID
        # space
        vg_args = ["vg", "ids", "-j"] + vg_names
        
        # It modifies the files in place
        options.drunner.call(job, [vg_args], work_dir=work_dir)
        
        cat_args = ["cat"] + vg_names
        
        with job.fileStore.writeGlobalFileStream() as (vg_handle, vg_id):
            # Make a merged vg with cat
            options.drunner.call(job, [cat_args], outfile=vg_handle,
                work_dir=work_dir)
            
        return vg_id
        
def shift_vgs_job(job, options, vg_ids):
    """
    Given a list of VG graph IDs, reassign the node IDs in the vg graphs so that
    they do not conflict. Returns a list of modified graph file IDs.
    
    TODO: use vg ids -j which is indeed in place.
    
    """
    
    if len(vg_ids) < 2:
        # Nothing to do because there are no possible conflicts with less than 2
        # graphs
        return vg_ids
    
    # First get all the ID pairs for all the graphs
    id_promises = [ToilPromise.wrap(job.addChildJobFn(toilvgfacade.id_range_job,
        options, vg)) for vg in vg_ids]
        
    # Then make jobs to shift all the graphs up. TODO: we could make this more
    # parallel by starting to shift the early graphs before we know how big the
    # late graphs are, with a sort of cascade.
    def shift_all(job, ranges):
        # This is an inline job to actually do all the shifting.
        # This holds the shift job RVs
        shift_results = []
        # This holds the max used ID
        used_ids = 0
        
        RealtimeLogger.info("Found ranges: {}".format(ranges))
        
        for (start, stop), vg_id in itertools.izip(ranges, vg_ids):
            # For each graph and its ID range
            
            # How much should it shift? Note that if we start at the max used ID
            # we need to shift by 1. Don't bother shifting down.
            shift = max(used_ids - start + 1, 0)
            
            # Shift it up
            shift_results.append(job.addChildJobFn(
                toilvgfacade.id_increment_job, options, vg_id, shift).rv())
            
            # What's the last ID we end up using?
            used_ids = stop + shift
            
        # Give back all the shifted graph IDs
        return shift_results
            
        
    # Return the result which is the list of shifted VG graph file IDs.
    return ToilPromise.all(id_promises).then_job_fn(shift_all).unwrap_result()
    
def explode_vg_job(job, options, plan, vg_id):
    """
    Explode the given vg file into connected component subgraphs.
    
    Returns a dict from subgraph ID to a list of path names present in the
    subgraph.
    """
    
    RealtimeLogger.info("Exploding graph {}".format(vg_id))
    
    # Define a work_dir for Docker
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the VG file
    vg_name = "graph.vg"
    job.fileStore.readGlobalFile(vg_id, os.path.join(work_dir, vg_name))
    
    # Set a directory for parts
    part_dir = "parts"
    
    # Reserve a report file
    report_filename = "report.tsv"
    
    # Set up args for the VG explode call
    vg_args = ["vg", "explode", vg_name, part_dir]
    
    # Run the call and capture the output
    options.drunner.call(job, [vg_args],
        outfile=open(os.path.join(work_dir, report_filename), "w"),
        work_dir=work_dir)
    
    RealtimeLogger.info("Exploded and put report in {}".format(
        os.path.join(work_dir, report_filename)))
    
    # We'll populate this with the paths that live in each part
    paths_by_part = {}
    
    for line in open(os.path.join(work_dir, report_filename)):
        if line == "":
            # Skip blank lines
            continue
        # Remove newlines and split each line on tabs.
        parts = line.strip().split("\t")
        
        # Validate the part
        options.drunner.call(job, [["vg", "validate", parts[0]]],
            work_dir=work_dir)
        
        # Save the file (first entry) to the filestore
        part_id = job.fileStore.writeGlobalFile(
            os.path.join(work_dir, parts[0]))
        
        # Remember that all of its paths belong to it. Make sure to trim tabs
        # and newlines.
        paths_by_part[part_id] = parts[1:]
        
        RealtimeLogger.info("Placed {} paths in part {}".format(
            len(parts) - 1, parts[0]))
        
    assert len(paths_by_part) > 0
        
    return paths_by_part
        
def add_variants_to_parts_job(job, options, plan, paths_by_part):
    """
    Given a dict from vg file ID to the paths in that VG file, work out the VCF
    files for each VG file and run jobs to add their variants. Return a dict of
    the same shape (vg graph ID to list of path names), but with the graph IDs
    replaced by IDs fro graphs with the variants added.
    
    """
    
    # This will hold VCF IDs by the path name they apply to
    vcfs_by_path = collections.defaultdict(list)
    
    for chromosome, vcf_id in plan.for_each_vcf_id_by_chromosome():
        # Build a dict from path name to relevant VCFs
        
        if not hasattr(options, 'add_chr') or options.add_chr:
            # Make sure to convert vcf names like "22" to UCSC-style "chr22".
            chromosome = "chr" + chromosome
            
        vcfs_by_path[chromosome].append(vcf_id)
    
    # Make a list of ToilPromises for updated vg graph file IDs
    updated_graphs = []
    
    # And a list of the lists of path names that go with them
    path_lists = []
    
    for vg_id, paths in paths_by_part.iteritems():
        # For each VG graph
        
        # Grab the applicable VCFs, flattening across all paths in the VG
        vcfs = list({vcf_id for path in paths for vcf_id in vcfs_by_path[path]})
        
        RealtimeLogger.info("VCFs: {} for paths: {}".format(vcfs, paths))
        
        # Add those VCFs to this VG
        child = job.addChildJobFn(add_variants_job, options, plan, vg_id,
            vcfs, cores=2, memory="100G", disk="50G")
            
        updated_graphs.append(ToilPromise.wrap(child))
        
        # And remember the paths in the graph
        path_lists.append(paths)
        
    # Wait for all the graphs to be done, use their IDs as dict keys in a dict
    # of path name lists, and return that.
    return ToilPromise.all(updated_graphs
        ).then(lambda graph_ids: dict(itertools.izip(graph_ids, path_lists))
        ).unwrap_result()
        
        
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
    merge_child = job.addChildJobFn(merge_vgs_job, options,
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
    
    # Define a work directory to mount in Dockers
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download all the VCFs
    vcf_filenames = []
    for i, vcf_id in enumerate(vcf_ids):
        # We need to download each VCF
        
        if plan.get_index_id(vcf_id) is not None:
            # This VCF is indexed, so it must be gzipped
            
            # Download the VCF under a .gz name
            vcf_filename = "vcf{}.vcf.gz".format(i)
            job.fileStore.readGlobalFile(vcf_id,
                os.path.join(work_dir, vcf_filename))
            
            # Download its index next to it
            job.fileStore.readGlobalFile(plan.get_index_id(vcf_id),
                os.path.join(work_dir, vcf_filename + ".tbi"))
            
            # Make sure the index is newer than the VCF
            os.utime(os.path.join(work_dir, vcf_filename + ".tbi"), None)
            
        else:
            # If it's not indexed, assume it's not compressed either.
            vcf_filename = "vcf{}.vcf".format(i)
            job.fileStore.readGlobalFile(vcf_id,
                os.path.join(work_dir, vcf_filename))
            
        # Save the filename
        vcf_filenames.append(vcf_filename)
        
        RealtimeLogger.info("Downloaded VCF {} as {}".format(vcf_id,
            os.path.join(work_dir, vcf_filename)))
    
    # Download the input graph
    vg_filename = "graph.vg"
    job.fileStore.readGlobalFile(vg_id, os.path.join(work_dir, vg_filename))
    
    # Validate it
    options.drunner.call(job, [["vg", "validate", vg_filename]],
        work_dir=work_dir)
    
    # Set up a command to add the VCFs to the graph
    vg_args = ["vg", "add", vg_filename]
    
    for vcf_filename in vcf_filenames:
        vg_args.append("-v")
        vg_args.append(vcf_filename)
        
        
    if options.add_chr:
        # Make sure to convert vcf names like "22" to UCSC-style "chr22".
        for base_name in plan.for_each_chromosome():
            vg_args.append("-n")
            vg_args.append("{}=chr{}".format(base_name, base_name))
    
    RealtimeLogger.info("Adding {} VCFs to vg graph...".format(
        len(vcf_filenames)))
    
    with job.fileStore.writeGlobalFileStream() as (vg_handle, new_vg_id):
        # Stream new graph to the filestore
        options.drunner.call(job, [vg_args], outfile=vg_handle,
            work_dir=work_dir)
    
    return new_vg_id
    
def hgvm_package_job(job, options, vg_id, xg_id, gcsa_pair, manifest_id):
    """
    Given a VG file ID, and XG file ID, a pair of GCSA and LCP file IDs, and a
    "serialized manifest ID, produce a Directory with "hgvm.vg", "hgvm.xg",
    ""hgvm.gcsa", hgvm.gcsa.lcp", and "hgvm.json" filled in.
    
    Runs as a Toil job so it will be able to unpack the pair.
    
    """
    
    # Package up the file IDs into a Directory
    hgvm = Directory()
    hgvm.add("hgvm.vg", vg_id)
    hgvm.add("hgvm.xg", xg_id)
    hgvm.add("hgvm.gcsa", gcsa_pair[0])
    hgvm.add("hgvm.gcsa.lcp", gcsa_pair[1])
    hgvm.add("hgvm.json", manifest_id)
    
    return hgvm
    
def graphs_to_hgvm_job(job, options, vg_ids, primary_paths=None):
    """
    Convert a list of vg graphs into an HGVM directory with its indexes.
    
    Can take a list of primary path names (or "") for each vg graph, for use
    with GCSA indexing. Those paths' nodes and edges are not pruned from the
    index for their graphs.
    
    Used for preparing a graph for realignment.
    
    Handles vg graph merging, and produces the hgvm manifest "hgvm.json".
    """
    
    # Write the manifest, which will include the primary paths.
    # We can add more facts here if we want.
    # TODO: hashes or something?
    manifest = Manifest({
        "hgvm_manifest_version": "0.1",
        "build_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "primary_paths": primary_paths,
        "uuid": str(uuid.uuid4())
    })
    
    # Merge all the graphs because we know that actually works
    merge_job = job.addChildJobFn(merge_vgs_job, options, vg_ids,
        cores=1, memory="100G", disk="20G")
    
    # Make the indexes using the individual graph files
    xg_job = merge_job.addFollowOnJobFn(toilvgfacade.xg_index_job, options,
        [merge_job.rv()])
    # Forward the primary paths on to the GCSA indexer, if specified, so they
    # can always be included.
    gcsa_job = merge_job.addFollowOnJobFn(toilvgfacade.gcsa_index_job, options,
        [merge_job.rv()], primary_paths)
    
    # Make a packaging job to package along with the single concatenated graph
    hgvm_job = xg_job.addFollowOnJobFn(hgvm_package_job, options,
        merge_job.rv(), xg_job.rv(), gcsa_job.rv(),
        manifest.save(job.fileStore))
    gcsa_job.addFollowOn(hgvm_job)
    merge_job.addFollowOn(hgvm_job)
    
    # Return its return value
    return hgvm_job.rv()
    
def hgvm_build_job(job, options, plan):
    """
    Toil job for building the HGVM. Returns a Directory object holding
    "hgvm.json", "hgvm.vg", "hgvm.xg", "hgvm.gcsa", and "hgvm.gcsa.lcp".
    
    """
    
    if plan.hgvm_directory is not None:
        # We have something to just import
        RealtimeLogger.info("Using imported HGVM")
        
        if options.dump_hgvm is not None:
            # Dump it again. TODO: unify dump logic
            plan.hgvm_directory.dump(job.fileStore, options.dump_hgvm)
        
        return plan.hgvm_directory
    
    # Otherwise we have to build it
    RealtimeLogger.info("Building HGVM")
    
    # Make a child to convert the HALs and merge the VGs
    vg_job = job.addChildJobFn(hals_and_vgs_to_vg_job, options, plan,
        list(plan.for_each_hal()), list(plan.for_each_base_vg()),
        cores=1, memory="2G", disk="1G")
        
    # Then split the VG up into connected components again
    explode_job = vg_job.addFollowOnJobFn(explode_vg_job, options, plan,
        vg_job.rv(),
        cores=1, memory="100G", disk="100G")
        
    # And add all the appropriate VCFs to each. Get a dict from updated graph
    # file ID to the list of paths it contains.
    add_job = explode_job.addFollowOnJobFn(add_variants_to_parts_job, options,
        plan, explode_job.rv(),
        cores=1, memory="4G", disk="4G")
        
    # Unpack that dict into keys and shortest path names (or empty strings)
    def reswizzle_path_names(paths_by_graph):
        # Return a pair of lists of graph IDs and corresponding inferred primary
        # path names. The primary path names are the shortest-named paths in
        # each graph, or "" if the graph has no paths.
        
        graph_ids = []
        primary_paths = []
        
        for graph_id, path_names in paths_by_graph.iteritems():
            # Order all the graphs
            graph_ids.append(graph_id)
            if len(path_names) == 0:
                # No primary path available for this chunk
                primary_paths.append("")
            else:
                # Grab the shortest name, sorting next by string value
                primary_paths.append(sorted(sorted(path_names), key=len)[0])
                
        RealtimeLogger.info("Guessed primary paths: {}".format(primary_paths))
        
        # Return the graphs and their guessed primary paths in the same order
        return graph_ids, primary_paths
    # Actually run that on the graph-id-to-name-list dict and get the promise
    reswizzle_promise = ToilPromise.wrap(add_job).then(reswizzle_path_names)
        
    # And then a job to merge and index and package it, using the primary paths
    # we just guessed for indexing. This job runs after the promise's job.
    hgvm_job = reswizzle_promise.addFollowOnJobFn(
        graphs_to_hgvm_job, options, reswizzle_promise.unwrap_result(0),
        primary_paths=reswizzle_promise.unwrap_result(1),
        cores=1, memory="4G", disk="4G")
    
    if options.dump_hgvm is not None:
        # Dump the hgvm directory from the filestore to local disk
        def dump(job, hgvm):
            hgvm.dump(job.fileStore, options.dump_hgvm)
        ToilPromise.wrap(hgvm_job).then_job_fn(dump)
    
    # Return the Directory with the packaged HGVM
    return hgvm_job.rv()
    
def concat_job(job, options, file_ids):
    """
    Merge zero or more VG protobuf files into one by concatenation. Returns the
    merged file ID.
    
    """
    
    with job.fileStore.writeGlobalFileStream() as (cat_handle, cat_id):
        # Make one merged file
        
        for part_id in file_ids:
            # For each part file
            with job.fileStore.readGlobalFileStream(part_id) as part_handle:
                # Open it
                
                # And stream it to the combined file
                shutil.copyfileobj(part_handle, cat_handle)
    
    return cat_id
    
def split_records_job(job, options, file_ids, lines_per_record=4,
    records_per_part=100000):
    """
    Split files that have a certain number of lines per record. Takes a list of
    file IDs, and a lines per record count. The default is set up for FASTQ
    files that have 4 lines per record.
    
    Each file in the list will be split at the same number of records, given by
    records_per_part.
    
    Returns a list of lists of corresponding parts of different input files. So
    if you pass in two files, the first parts of the files are in the first
    list, the second parts of the files are in the second list, and so on.
    
    TODO: replace with a real FASTQ parser that can handle valid FASTQs that
    wrap over more than 4 lines.
    """
    
    if len(file_ids) == 0:
        # Super easy to split no files.
        return []
    
    # Define a work_dir to work in for Docker
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the files
    file_names = ["file{}".format(i) for i in xrange(len(file_ids))]
    for file_id, file_name in itertools.izip(file_ids, file_names):
        # Each one goes undar a name we assigned in work_dir
        job.fileStore.readGlobalFile(file_id, os.path.join(work_dir, file_name))
    
    # Measure the first file
    line_count = int(options.drunner.call(job, [["wc", "-l", file_names[0]],
        ["cut", "-f", "1", "-d", " "]],
        check_output=True, work_dir=work_dir))
    
    # There may not be any partial records
    assert(line_count % lines_per_record == 0)
    
    # Work out how many actual parts we need
    whole_parts, partial_parts = divmod(line_count,
        lines_per_record * records_per_part)
    part_count = whole_parts + int(bool(partial_parts))
    
    RealtimeLogger.info("Splitting {} into {} parts of {} lines each...".format(
        file_ids, part_count, lines_per_record * records_per_part))
    
    # This holds a list for each chunk number of all the corresponding chunks
    # from all files, as file IDs
    part_lists = []
    
    # This holds the file IDs for the chunk currently being worked on
    part_list = []
    
    # Open all the files to chunk
    input_handles = [open(os.path.join(work_dir, x)) for x in file_names]
    
    for i in xrange(part_count):
        # For each part we need to do
    
        for input_handle in input_handles:
            # For each file that's open
            with job.fileStore.writeGlobalFileStream() as (chunk_handle,
                chunk_id):
                # Open a chunk file to fill in
                for line in itertools.islice(input_handle,
                    records_per_part * lines_per_record):
                    # Read to the end or to the requested number of records
                    
                    # And put them in the chunk file
                    chunk_handle.write(line)
                
                # Finish the chunk and stick it in its list
                part_list.append(chunk_id)
        
        # Now we have a whole set of corresponding parts, so send it out
        RealtimeLogger.info("Part {}: {}".format(i, part_list))
        part_lists.append(part_list)
        part_list = []
    
    for input_handle in input_handles:
        # Make sure all the files are empty now. They all should have had the same length.
        assert(not input_handle.read(1))
        
    # Return the list of lists of corresponding chunks
    return part_lists
    
def align_to_hgvm_chunked_job(job, options, hgvm, fastqs=[],
    sequences=None):
    """
    Align either a single or a pair of paired FASTQs, or a file of sequences,
    one per line, to an HGVM represented as a Directory.
    
    Only one of sequences and fastqs will be used.
    
    Return the GAM file ID.
    
    Aligns in multiple chunks
    """
    
    if sequences is not None:
        # We're using sequence files
        
        split_job = job.addChildJobFn(split_records_job, options, [sequences],
            lines_per_record=1, cores=1, memory="4G", disk="100G")
        
    else:
        # We must be using fastqs
        split_job = job.addChildJobFn(split_records_job, options, fastqs,
            lines_per_record=4, cores=1, memory="4G", disk="100G")
          
    # We define this little inline job to actually do the alignments once the
    # splitting happens
    def do_alignments(job, parts):
        # This holds the .rv()s for the aligned GAM files
        rvs = []
        for part_list in parts:
            if sequences is not None:
                # Use the one-element parts lists as sequence files
                rvs.append(job.addChildJobFn(align_to_hgvm_job, options, hgvm,
                    sequences=part_list[0],
                    cores=32, memory="200G", disk="200G").rv())
            else:
                # Use the multi-element parts lists as FASTQs
                rvs.append(job.addChildJobFn(align_to_hgvm_job, options, hgvm,
                    fastqs=part_list,
                    cores=32, memory="200G", disk="200G").rv())
        # Return all the aligned parts
        return rvs
            
    # After splitting, do all the alignments, then merge the GAMs and return
    # that
    return ToilPromise.wrap(split_job).then_job_fn(do_alignments).then_job_fn(
        concat_job, options).unwrap_result()
    
    
    
    
def align_to_hgvm_job(job, options, hgvm, fastqs=[], sequences=None):
    """
    Align either a single or a pair of paired FASTQs, or a file of sequences,
    one per line, to an HGVM represented as a Directory.
    
    Return the GAM file ID.
    
    Does the whole alignment as a single job.
    """
    
    # Make a work directory for Docker
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the HGVM
    hgvm_dir = "hgvm"
    hgvm.download(job.fileStore, os.path.join(work_dir, hgvm_dir))
    
    # Prepare some args to align to the HGVM
    vg_args = ["vg", "map", "-t", "32", "-x", os.path.join(hgvm_dir, "hgvm.xg"),
        "-g", os.path.join(hgvm_dir, "hgvm.gcsa")]
        
    for i, fastq_id in enumerate(fastqs):
        # Download and point at all the FASTQs
        fastq_filename = "fastq{}.fastq".format(i)
        job.fileStore.readGlobalFile(fastq_id,
            os.path.join(work_dir, fastq_filename))
        vg_args.append("--fastq")
        vg_args.append(fastq_filename)
        
    if sequences is not None:
        # Download and point to the sequences file
        sequences_filename = "sequences.txt"
        job.fileStore.readGlobalFile(sequences,
            os.path.join(work_dir, sequences_filename))
        vg_args.append("--reads")
        vg_args.append(sequences_filename)
    
    # TODO: configure multimapping and stuff
    
    RealtimeLogger.info("Aligning FASTQs: {} and sequences: {} to HGVM".format(
        fastqs, sequences))
    
    with job.fileStore.writeGlobalFileStream() as (gam_handle, gam_id):
        # Align and stream GAM to the filestore
        options.drunner.call(job, [vg_args], outfile=gam_handle,
            work_dir=work_dir)
    
    return gam_id
    
def pileup_on_hgvm_job(job, options, hgvm, gam_id):
    """
    Pile up the given GAM file on the given packaged HGVM. Returns the pileup
    file ID.
    
    """

    # Make a work directory for Docker
    work_dir = job.fileStore.getLocalTempDir()

    # Download just the vg
    vg_filename = "graph.vg"
    job.fileStore.readGlobalFile(hgvm.get("hgvm.vg"),
        os.path.join(work_dir, vg_filename))
    
    # Download the GAM
    gam_filename = "aligned.gam"
    job.fileStore.readGlobalFile(gam_id, os.path.join(work_dir, gam_filename))
    
    # Set up the VG run
    vg_args = ["vg", "pileup", vg_filename, gam_filename]
    
    # TODO: configure filtering and stuff
    
    try:
        with job.fileStore.writeGlobalFileStream() as (pileup_handle, pileup_id):
            # Pile up and stream the pileup to the file store
            options.drunner.call(job, [vg_args], outfile=pileup_handle,
                work_dir=work_dir)
    except Exception as e:
        logging.error("Failed. Dumping files.")
        
        hgvm.add("aligned.gam", gam_id)
        hgvm.dump(job.fileStore, "/home/anovak/dump/{}".format(
            random.randint(0, 1000000000)))
        raise e
    
    
    return pileup_id
    
def call_on_hgvm_job(job, options, hgvm, pileup_id, vcf=False):
    """
    Given a packaged HGVM Directory and a pileup file ID, produce variant calls
    in Locus format. Returns a Directory with the Locus data as "calls.loci" or
    the VCF data as "calls.vcf", and the augmented graph as "augmented.vg".
    
    """


    if hgvm.has_file("hgvm.json"):
        # Load the manifest, which lists the primary paths, which we want to
        # know for calling
        manifest = Manifest.load(job.fileStore, hgvm.get("hgvm.json"))
    else:
        # Make a fake Manifest with no primary paths. vg call had better be able
        # to guess them itself.
        RealtimeLogger.warn("Creating fake manifest with no primary paths!")
        manifest = Manifest()
        
    # Kick off a job to do the actual calling, and return the vcf/loci and the
    # augmented graph.
    call_job = job.addChildJobFn(toilvgfacade.vg_call_job, options,
        hgvm.get("hgvm.vg"), pileup_id, vcf=vcf,
        primary_paths=manifest.get("primary_paths", []))
        
    # Package the results into a directory and return it.
    return ToilPromise.wrap(call_job
        ).then(lambda (output_id, augmented_id): Directory({
            "calls.vcf" if vcf else "calls.loci": output_id,
            "augmented.vg": augmented_id
        })).unwrap_result()
    
def subset_graph_job(job, options, eval_plan, call_directory):
    """
    Given a Directory with an "augmented.vg" and a "calls.loci", subset the
    augmented graph to the sample graph defined by the Loci, and return the new
    sample graph's file ID.
    """
    
    # Define a work_dir so Docker can work
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the calling results
    call_dir = "call"
    call_directory.download(job.fileStore, os.path.join(work_dir, call_dir))
    
    # Set up the VG run
    vg_args = ["vg", "mod", "--sample-graph", call_dir + "/calls.loci",
        call_dir + "/augmented.vg"]
        
    with job.fileStore.writeGlobalFileStream() as (sample_handle, sample_id):
        # Stream the sample graph to the file store
        options.drunner.call(job, [vg_args], outfile=sample_handle,
            work_dir=work_dir)
        
    return sample_id
    
def realignment_stats_job(job, options, eval_plan, vg_id, gam_id):
    """
    Use vg stats to get a report for the given alignments to the given graph.
    
    """
    
    # Define a work_dir so Docker can work
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the VG file
    vg_name = "graph.vg"
    job.fileStore.readGlobalFile(vg_id, os.path.join(work_dir, vg_name))
    
    # Download the GAM file
    gam_name = "aligned.gam"
    job.fileStore.readGlobalFile(gam_id, os.path.join(work_dir, gam_name))
    
    # Set up the VG run
    vg_args = ["vg", "stats", "--verbose", "--alignments", gam_name, vg_name]
        
    with job.fileStore.writeGlobalFileStream() as (stats_handle, stats_id):
        # Stream the stats to the file store
        options.drunner.call(job, [vg_args], outfile=stats_handle,
            work_dir=work_dir)
        
    return stats_id
    
def parse_realignment_stats_job(job, options, eval_plan, stats_id):
    """
    Parse the output of vg stats on a GAM. Returns two dicts.
    
    The first is a dict from stat name ("Insertions", "Deletions",
    "Substitutions", "Softclips") to (bp, event count) pairs.
    
    The second is a dict from node status name ("Unvisited nodes", "Single-
    visited nodes") to (bp, node count) pairs.
    
    """
    
    # Make the dict to fill in with stats
    stats = {}
    
    # Make another dict to fill in with node visit status
    node_status = {}
    
    # We have some regexes to catch lines we care about
    
    # These are alignment events
    stats_regex = re.compile(R"(.+): (\d+) bp in (\d+) read events")
    
    # These are node statuses
    node_regex = re.compile(R"(.+): (\d+)/(\d+) \((\d+) bp\)")
    
    with job.fileStore.readGlobalFileStream(stats_id) as stats_file:
        for line in stats_file:
            # For every line, see if it's a normal stat-giving line
            stats_match = stats_regex.match(line.strip())
            if stats_match:
                # Put it in the dict as (bp, event count)
                stats[stats_match.group(1)] = (int(stats_match.group(2)),
                    int(stats_match.group(3)))
                    
            # Also see if it describes node visit statsu
            node_match = node_regex.match(line.strip())
            if node_match:
                # Put it in the dict as (bp, node count), which is reverse from
                # how it comes in.
                node_status[node_match.group(1)] = (int(node_match.group(3)),
                    int(node_match.group(2)))
            
    # Return the completed dictionaries
    return stats, node_status
            
def measure_graph_job(job, options, vg_id):
    """
    Measue the given graph with vg stats. Return a pair of (bp size, node
    count).
    """
    
    # Define a work_dir so Docker can work
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download just the vg
    vg_filename = "graph.vg"
    job.fileStore.readGlobalFile(vg_id, os.path.join(work_dir, vg_filename))
    
    # Set up the VG run
    vg_args = ["vg", "stats", "--size", "--length", vg_filename]
    
    # Call and stream the Locus data to the file store
    stats = options.drunner.call(job, [vg_args], check_output=True,
        work_dir=work_dir)
    
    # Fill in these stat variables
    bp_length = None
    node_count = None
    
    for line in stats.split("\n"):
        # Break each line on the tab
        parts = line.split("\t")
        
        if len(parts) < 2:
            # Skip tiny nonsense lines
            continue
        
        if parts[0] == "length":
            # Parse the total length
            bp_length = int(parts[1])
            
        if parts[0] == "nodes":
            # Parse the total node count
            node_count = int(parts[1])
            
    # Return the two metrics of graph size
    return bp_length, node_count
    
def stats_tsv_job(job, options, eval_plan, stats):
    """
    Takes some statistics, in the form of a pair.
    
    The first element is a dict from condition to pair of dicts, one of edit
    operation stats and one of node status stats. Each dict entry is a pair of
    (bp, event count).
    
    The second element is a pair (bp length, node count) for the control graph,
    used for normalizing all the stats in terms of events or bases per primary
    reference base. It may be None, in which case we should not normalize.
    
    Returns a Directory of TSV files, named after edit operation stats.
    """
    
    # Unpack the argument
    (by_condition, control_size) = stats
    
    # Define a directory to fill
    directory = Directory()
    
    for edit in ("Insertions", "Deletions", "Substitutions", "Softclips"):
        # For each kind of edit operation, we'll make a TSV in the file store
    
        with job.fileStore.writeGlobalFileStream() as (tsv_handle, tsv_id):
            writer = tsv.TsvWriter(tsv_handle)
        
            for condition, (edit_stats, node_stats) in by_condition.iteritems():
                # For each condition
                
                # Grab the edited bp count
                edited = edit_stats[edit][0]
        
                if control_size is not None:
                    # Normalize to fraction of primary path bases, if we're
                    # normalizing
                    edited /= float(control_size[0])
                
                # Add a TSV record for this condition
                writer.line(condition, edited)
                
            # Then put the TSV in the directory after we put in all the
            # conditions
            directory.add("{}.tsv".format(edit), tsv_id)
            
    return directory
        
    
def hgvm_eval_job(job, options, eval_plan, hgvm):
    """
    Toil job for evaluating the HGVM. Takes a packaged HGVM Directory, and a
    plan for evaluating it.
    
    Returns the evaluation results as a Directory.
    
    """
    
    RealtimeLogger.info("Evaluating HGVM")
    
    if len(eval_plan.get_fastq_ids()) > 0:
        # If we have reads, we can do an evaluation
        
        # We'll fill this in with all the graphs we want to test, as packaged
        # HGVM directories or .rv() promises for them.
        graphs_to_test = {"experimental": hgvm}
        
        # This will hold promises for evaluation directories, by condition name.
        eval_promises = {}
        
        # This will hold promises for parsed stat dict pairs from
        # parse_realignment_stats_job, by condition name.
        stats_promises = {}
        
        # This holds a ToilPromise for either the control graph's length in bp,
        # or None if no control graph is specified.
        control_length_promise = None
        
        control_hgvm = eval_plan.get_packaged_control()
        
        if control_hgvm is not None:
            # Add in the control. TODO: assumes it is already indexed. We should
            # check and index it if it's not.
            graphs_to_test["control"] = control_hgvm
            
            # Get a promise for the control's length
            control_length_promise=ToilPromise.wrap(job.addChildJobFn(
                measure_graph_job, options,
                control_hgvm.get("hgvm.vg"),
                cores=1, memory="50G", disk="50G"))
        else:
            # Say we can't know the control's length
            control_length_promise = ToilPromise.resolve(job, None)
        
        for condition, package in graphs_to_test.iteritems():
            # Looping over the experimental and the control...
        
            align_job = job.addChildJobFn(align_to_hgvm_chunked_job, options,
                package, fastqs=eval_plan.get_fastq_ids(),
                cores=16, memory="50G", disk="100G")
            # Put it in a Directory in a promise
            align_promise = ToilPromise.wrap(align_job).then(
                lambda id: Directory({"aligned.gam": id}))
            
            # Then do the pileup
            pileup_job = align_job.addFollowOnJobFn(pileup_on_hgvm_job, options,
                package, align_job.rv(),
                cores=1, memory="50G", disk="100G")
            # Put it in a Directory in a promise
            pileup_promise = ToilPromise.wrap(pileup_job).then(
                lambda id: Directory({"pileup.vgpu": id}))
            
            # Then do the calling
            call_job = pileup_job.addFollowOnJobFn(call_on_hgvm_job, options,
                package, pileup_job.rv(),
                cores=1, memory="100G", disk="50G")
            # And a promise for the call results directory.
            call_promise = ToilPromise.wrap(call_job)
            
            # Then subset the graph and get the new VG id and put it in a
            # directory. Then do the calling
            subset_job = call_job.addFollowOnJobFn(subset_graph_job, options,
                eval_plan, call_job.rv(),
                cores=1, memory="50G", disk="50G")
            # And a promise for the call results directory.
            subset_promise = ToilPromise.wrap(subset_job).then(
                lambda id: Directory({"sample.vg": id}))
                
            # Make a list of all the output directory promises to squash together.
            dir_promises = [align_promise, pileup_promise, call_promise,
                subset_promise]
                
            if eval_plan.get_eval_sequences_id() is not None:
                # We should realign these sequences
                
                # First reindex and package the sample itself as an HGVM
                package_job = subset_job.addFollowOnJobFn(
                    graphs_to_hgvm_job, options, [subset_job.rv()])
                
                # Then realign and grab the realigned GAM
                realign_job = package_job.addFollowOnJobFn(
                    align_to_hgvm_chunked_job, options, package_job.rv(),
                    sequences=eval_plan.get_eval_sequences_id(),
                    cores=16, memory="50G", disk="100G")
                # And put that in a directory
                realign_promise = ToilPromise.wrap(realign_job).then(
                    lambda id: Directory({"realigned.gam": id}))
                # And merge it with the others
                dir_promises.append(realign_promise)
                    
                # Then get the stats
                stats_job = realign_job.addFollowOnJobFn(realignment_stats_job,
                    options, eval_plan, subset_job.rv(), realign_job.rv(),
                    cores=1, memory="50G", disk="100G")
                # And put them in a directory
                stats_promise = ToilPromise.wrap(stats_job).then(
                    lambda id: Directory({"stats.tsv": id}))
                # And merge it with the others
                dir_promises.append(stats_promise)
                
                # Then parse them
                stats_parse_job = stats_job.addFollowOnJobFn(
                    parse_realignment_stats_job,
                    options, eval_plan, stats_job.rv())
                # And make a promise and save it
                stats_promises[condition] = ToilPromise.wrap(stats_parse_job)

            # Collect the result directories and merge them. Store a promise for
            # that directory under the condition name we are doing.
            eval_promises[condition] = ToilPromise.all(dir_promises).then(
                lambda dirs: dirs[0].merge(*dirs[1:]))
                
        if (eval_plan.get_eval_sequences_id() is not None and
            graphs_to_test.has_key("control")):
            # Also consider realigning just against the control graph, without
            # calling
            condition = "allref"
            
            # Control is already indexed and packaged.
            
            # Realign and grab the realigned GAM
            realign_job = job.addChildJobFn(align_to_hgvm_chunked_job,
                options, graphs_to_test["control"],
                sequences=eval_plan.get_eval_sequences_id(),
                cores=16, memory="50G", disk="100G")
                
            # Then get the stats
            stats_job = realign_job.addFollowOnJobFn(realignment_stats_job,
                options, eval_plan, graphs_to_test["control"].get("hgvm.vg"),
                realign_job.rv(),
                cores=1, memory="50G", disk="100G")
                
            # Put both in a directory
            realign_dir_promise = ToilPromise.all({
                "realigned.gam": ToilPromise.wrap(realign_job),
                "stats.tsv": ToilPromise.wrap(stats_job)
            }).then(lambda x: Directory(x))
                
            # And make a directory with just the stats for this condition
            eval_promises[condition] = realign_dir_promise
            
            # Then parse them
            stats_parse_job = stats_job.addFollowOnJobFn(
                parse_realignment_stats_job,
                options, eval_plan, stats_job.rv())
            # And make a promise and save it
            stats_promises[condition] = ToilPromise.wrap(stats_parse_job)
                
        # Now we have evaluated all the conditions. Make a plot TSV.
        
        # Gather all the stats dict pairs by condition, and the control length
        # (if any), and compute a Directory of TSVs using the appropriate Toil
        # job function.
        tsv_dir_promise = ToilPromise.all([ToilPromise.all(stats_promises),
            control_length_promise]).then_job_fn(stats_tsv_job,
                options, eval_plan)
        
        # Mount each condition's own directory in the eval directory.
        def mount_all(dirs_by_condition):
            root = Directory()
            for condition, directory in dirs_by_condition.iteritems():
                root.mount(condition, directory)
            return root
        eval_dir_promise = ToilPromise.all(eval_promises).then(mount_all)
        
        # Then stick in all the comparative TSVs, and return the result
        return ToilPromise.all([eval_dir_promise, tsv_dir_promise]).then(
            lambda args: args[0].mount("tsv", args[1])).unwrap_result()

    # Otherwise return this empty directory
    return Directory()
    
def measure_sv_recall_job(job, options, call_dir, truth_vcfs, sample_name,
    threshold=25, radius=25):
    """
    Given a directory with a "calls.vcf" in it, and a bunch of file IDs for
    truth VCFs that may contain some structural variants, calculate the portion
    of SVs in the truth set that have SVs called sufficiently close in the call
    VCF.
    
    threshold gives the minimum length change for a called variant to be
    structural, and radius gives how close a called SV has to be to an expected
    SV for the expected SV to be counted as recalled. Truth variants are
    identified as structural by the presence of an SVTYPE.
    
    Returns a Directory with a "recall.tsv" containing the statistics.
    """
    
    RealtimeLogger.info("Compare calls against {} truth VCFs".format(
        len(truth_vcfs)))
    
    # First load all the SV call positions and cover a range around the start in
    # an interval tree per contig
    contig_trees = collections.defaultdict(intervaltree.IntervalTree)
    
    with job.fileStore.readGlobalFileStream(call_dir.get("calls.vcf")) as calls:
        # Stream the file
        for record in vcf.Reader(calls):
            # For each call
            
            # Filters get parsed into a list, or None if it's "."
            if record.FILTER is not None and len(record.FILTER) > 0:
                # Skip filtered variants
                continue
            
            # Find out its max SVLEN length change
            svlen = max([abs(int(x)) for x in record.INFO.get("SVLEN", [0])])
            
            if svlen < threshold:
                # Too small to be an SV
                continue
                
            if record.samples[0]["GT"] in {"0/0", "0|0", "0"}:
                # Skip SVs that aren't actually called in the sample
                continue
            
            # The calls that come out already have "chr" if necessary.
            
            # Add an interval around the SV call
            contig_trees[record.CHROM].addi(record.POS - radius,
                record.POS + radius, record.POS)
                
    # Set up some stats to fill in
    found_svs = 0
    total_svs = 0
    
    # This maps from recalled position to positions that the recalling calls are
    # at
    recalling_positions = collections.defaultdict(set)
                
    # This holds a set of seen CHROM, POS pairs, so we don't count the duplicate
    # multiple-lengths-of-a-single-SV variants as the same variant.
    seen_truth_starts = set()
                
    # Then loop over the truth VCFs
    for vcf_id in truth_vcfs:
        with job.fileStore.readGlobalFileStream(vcf_id) as truth:
            reader = vcf.Reader(TransparentUnzip(truth))
            for record in reader:
                # For every truth VCF record
                
                if not record.INFO.has_key("SVTYPE"):
                    # Skip non-structural variants
                    continue
                    
                if record.FILTER is not None and len(record.FILTER) > 0:
                    # Skip filtered variants
                    continue
                    
                if record.genotype(sample_name)["GT"] in {"0/0", "0|0", "0"}:
                    # Skip SVs that aren't actually supposed to be in the sample
                    continue
                    
                # Has this been seen before?
                pos_key = (record.CHROM, record.POS)
                if pos_key in seen_truth_starts:
                    # This has been seen before (or another variant at the same
                    # site, which we're going to consider the same for our
                    # purposes).
                    continue
                seen_truth_starts.add(pos_key)
                    
                # This is an SV
                total_svs += 1
                
                # See what contig it's on
                chrom_name = record.CHROM
                if options.add_chr:
                    chrom_name = "chr" + chrom_name
                    
                if len(contig_trees[chrom_name][record.POS]) != 0:
                    # We hit near enough to a called variant, so this SV is
                    # recalled.
                    found_svs += 1
                    
                    for recalling in contig_trees[chrom_name][record.POS]:
                        recalling_positions[record.POS].add(recalling.data)
    
    # Compute the portion recalled and output a TSV with it and the totals
    portion_recalled = float(found_svs) / total_svs if total_svs > 0 else 0
    
    # Drop the TSV in a Directory and return it
    with job.fileStore.writeGlobalFileStream() as (tsv_handle, tsv_id):
        writer = tsv.TsvWriter(tsv_handle)
        writer.line("Found", found_svs)
        writer.line("Total", total_svs)
        writer.line("Recall", portion_recalled)
        
        for pos, recalled_by in recalling_positions.iteritems():
            # Dump all the recall relationships
            writer.comment("{}: {}".format(pos, ",".join(
                [str(x) for x in recalled_by])))
        
        return Directory({"recall.tsv": tsv_id})
    
def hgvm_sv_recall_job(job, options, plan, recall_plan, hgvm):
    """
    Toil job for evaluating an HGVM's structural variant calling recall rate.
    
    """
    
    # Make sure we have everything we would need to actually do this job
    
    if len(recall_plan.get_sample_names()) == 0:
        # Don't do anything!
        return Directory()
    
    # Otherwise, we fill this with promises for directories, one per sample
    sample_promises = []
    
    for sample_name in recall_plan.get_sample_names():
        # For each sample we want to compute recall for
    
        # Within each sample we promise Directories with evaluation files which
        # we are going to merge.
        dir_promises = []
    
        if recall_plan.get_vcf_id(sample_name) is not None:
            # We already have a VCF for this sample. Use it.
            
            # Just emit a Directory with those as the calls
            call_promise = ToilPromise.resolve(job,
                Directory({"calls.vcf": recall_plan.get_vcf_id(sample_name)}))
            dir_promises.append(call_promise)
        else:
            if recall_plan.get_gam_id(sample_name) is not None:
                # We can just use this GAM
                align_promise = ToilPromise.resolve(job, Directory(
                    {"aligned.gam": recall_plan.get_gam_id(sample_name)}))
                dir_promises.append(align_promise)
                # Do the pileup as a child
                pileup_job = job.addChildJobFn(pileup_on_hgvm_job, options,
                    hgvm, recall_plan.get_gam_id(sample_name),
                    cores=1, memory="50G", disk="100G")
            else:
                # We have to actually align some FASTQs    
                fastqs=recall_plan.get_fastq_ids(sample_name)
                if len(fastqs) == 0:
                    # We have a sample with no FASTQ and no GAM
                    raise RuntimeError("No reads available for recall evaluation "
                        "on {}".format(sample_name))
                
                # Align the reads
                align_job = job.addChildJobFn(align_to_hgvm_chunked_job, options,
                    hgvm, fastqs=fastqs,
                    cores=16, memory="50G", disk="100G")
                # Put it in a Directory in a promise
                align_promise = ToilPromise.wrap(align_job).then(
                    lambda id: Directory({"aligned.gam": id}))
                dir_promises.append(align_promise)
                    
                # Then do the pileup as a follow-on
                pileup_job = align_job.addFollowOnJobFn(pileup_on_hgvm_job, options,
                    hgvm, align_job.rv(),
                    cores=1, memory="50G", disk="100G")
            # Put it in a Directory in a promise
            pileup_promise = ToilPromise.wrap(pileup_job).then(
                lambda id: Directory({"pileup.vgpu": id}))
            dir_promises.append(pileup_promise)
        
            # Then do the calling
            call_job = pileup_job.addFollowOnJobFn(call_on_hgvm_job, options,
                hgvm, pileup_job.rv(), vcf=True,
                cores=1, memory="100G", disk="50G")
            # And a promise for the call results directory.
            call_promise = ToilPromise.wrap(call_job)
            dir_promises.append(call_promise)
        
        # Collect the expected variants from all the input VCFs for the
        # chromosomes we did
        truth_vcfs = []
        for chromosome in plan.for_each_chromosome():
            # TODO: limit to just the contigs we are actually doing, somehow...
            for vcf_id in plan.get_vcf_ids(chromosome):
                # Aggregate into a list of VCFs
                truth_vcfs.append(vcf_id)
        
        # Compare the expected variants to the observed variants and compute the
        # fraction recalled
        compare_job = call_promise.addFollowOnJobFn(measure_sv_recall_job,
            options, call_promise.unwrap_result(), truth_vcfs, sample_name,
            cores=1, memory="8G", disk="100G")
        # And add a promise for its results Directory
        compare_promise = ToilPromise.wrap(compare_job)
        dir_promises.append(compare_promise)
    
        # Promise all that in a Directory
        sample_promises.append(ToilPromise.all(dir_promises).then(
            lambda dirs: dirs[0].merge(*dirs[1:])))
         
    # Mount all the per-sample directories under their sample names   
    return ToilPromise.all(sample_promises
        ).then(lambda dirs: Directory().mount_all(
        dict(itertools.izip(recall_plan.get_sample_names(), dirs)))
        ).unwrap_result()
    
def main_job(job, options, plan, eval_plan, recall_plan):
    """
    Root job of the Toil workflow. Build and then evaluate an HGVM.
    
    Returns a fully-structured output Directory. All the master has to do is
    export it.
    
    """
    
    RealtimeLogger.info("Main job starting")
    
    # Make sure VG is available
    options.drunner.call(job, [["vg", "version"]])
    
    # And samtools
    # TODO: use vg to index FASTAs once, somehow.
    options.drunner.call(job, [["samtools", "--version"]])
    
    # And hal2vg (which can't be made to succeed with no real arguments...)
    #options.drunner.call(job, [["which", "hal2vg"]])
    # TODO: Dockerize hal2vg
    
    # Build the HGVM
    build_job = job.addChildJobFn(hgvm_build_job, options, plan,
        cores=1, memory="2G", disk="1G")
    
    # Then evaluate it with the realignment evaluation
    eval_job = build_job.addFollowOnJobFn(hgvm_eval_job, options, eval_plan,
        build_job.rv(),
        cores=1, memory="2G", disk="1G")
        
    # And with the SV recall evaluation
    sv_eval_job = build_job.addFollowOnJobFn(hgvm_sv_recall_job, options, plan,
        recall_plan, build_job.rv(),
        cores=1, memory="2G", disk="1G")
        
    # Stick all those directories in a folder structure and return that
    return ToilPromise.all({
        "hgvm": ToilPromise.wrap(build_job),
        "eval/realignment": ToilPromise.wrap(eval_job),
        "eval/sv": ToilPromise.wrap(sv_eval_job)
    }).then(lambda dirs: Directory().mount_all(dirs)).unwrap_result()
    
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
    
    # Add the drunner to the options and initialize stuff from the Toil VG
    # config
    toilvgfacade.initialize(options)
    
    # Start up Toil
    with Toil(options) as toil_instance:
        
        if toil_instance.options.restart:
            # We're re-running. Grab the root job return value from restart
            directory = toil_instance.restart()
        else:
            # Run from the top
        
            # Build the plan on the head node
            plan = create_plan(options.assembly_url, options.vcfs_url,
                options.vcf_url, options.vcf_contig, options.hal_url,
                options.base_vg_url, options.hgvm_url)
                
            # Also build the realignment evaluation plan on the head node
            eval_plan = create_eval_plan(options.control_graph_hgvm_url, 
                options.control_graph_url, options.control_graph_xg_url,
                options.control_graph_gcsa_url, options.sample_fastq_url,
                options.eval_sequences_url)
                
            # And the SV recall plan
            recall_plan = create_recall_plan(options.sv_sample_fastq_url,
                options.sv_sample_gam_url, options.sv_sample_vcf_url,
                options.sv_sample_name)
        
            # Import all the files from the plans. Now the plans will hold Toil
            # IDs for data files, and actual info for metadata files.
            plan.bake(lambda url: toil_instance.importFile(url))
            eval_plan.bake(lambda url: toil_instance.importFile(url))
            recall_plan.bake(lambda url: toil_instance.importFile(url))
        
            # Make a root job
            root_job = Job.wrapJobFn(main_job, options, plan, eval_plan,
                recall_plan,
                cores=1, memory="1G", disk="1G")
            
            # Run the root job and get the final output directory
            directory = toil_instance.start(root_job)
        
        # Export the results
        directory.export_to(
            lambda id, url: toil_instance.exportFile(id, url),
            options.out_url)
        
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
        
        
        
        
        
        
        
        
        
        

