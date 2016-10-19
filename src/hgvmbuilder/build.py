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

import toil
from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from .plan import ReferencePlan
from . import grcparser
from . import thousandgenomesparser
from .vcfrewriter import VcfRewriter

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
    
def rewrite_vcf_job(job, options, plan, vcf_id):
    """
    Given the ID of a vcf file (as well as the options and plan), return the ID
    of a rewritten VCF where chromosome names have been translated into
    accession.version format.
    """
    
    # Make a rewriter to rewrite the VCF
    rewriter = VcfRewriter(plan.get_name_translation())
    
    
    with job.fileStore.readGlobalFileStream(vcf_id) as input_stream:
        # Read the input file
        
        with job.fileStore.writeGlobalFileStream(cleanup=True) as \
            (output_stream, new_id):
            # Fill the output file
            
            # Rewrite VCF from one stream to the other
            rewriter.rewrite_stream(input_stream, output_stream)
        
    RealTimeLogger.info("Rewrote VCF to {}".format(new_id))
        
    # Return the ID for the rewritten, uncompressed, unindexed VCF
    return new_id
    
    # TODO: compress VCF with bgzip, index with tabix
    
def main_job(job, options, plan):
    """
    Root Toil job. Execute the plan.
    """

    # Farm out jobs to break up all the FASTAs by pertinent chromosome
    
    # Collect together file IDs for all the chromosomes for primary (where there
    # should be one) and the alts
    
    for chromosome_name, vcf_id in plan.for_each_vcf_id_by_chromosome():
        # In parallel we can take the VCFs and translate them to
        # accession.version. Then we have to re-compress and re-index.
        
        new_id = job.addChildJobFn(rewrite_vcf_job, options, plan, vcf_id).rv()
        
        # TODO: pass to re-indexing child.
        
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
            # We're re-running
            toil_instance.restart()
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
            toil_instance.start(root_job)
        
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
        
        
        
        
        
        
        
        
        
        

