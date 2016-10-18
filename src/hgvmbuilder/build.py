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

from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from .plan import ReferencePlan
from .ftputil import FTPOrFilesystemConnection

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
    parser.add_argument("--assembly_structure",
        default=("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCA_000001405.24_GRCh38.p9/"
        "GCA_000001405.24_GRCh38.p9_assembly_structure"),
        help="root of input assembly structure in GRC format")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def create_plan(assembly_structure):
    """
    Given an FTP or file url to the root of a GRC-format assembly_structure
    directory tree, produce a ReferencePlan describing that assembly.
    
    Recursively traverses the directory structure and finds the various
    FASTAs and metadata files and adds their URLs to the plan.
    """
    
    # Make the plan
    plan = ReferencePlan()
    
    # Open the URL
    connection = FTPOrFilesystemConnection(assembly_structure)
    
    for assembly_unit_root in connection.list_children(""):
        # For each assembly unit
        Logger.debug("Assembly unit: {}".format(assembly_unit_root))
        for scaffold_category in connection.list_children(assembly_unit_root):
            # For everything under that
            
            # Get the base name (last path component)
            basename = os.path.basename(scaffold_category)
            
            Logger.debug("Scaffold category: {}".format(basename))
            
            if basename == "placed_scaffolds":
                # These are redundant with the actual assembled contigs
                continue
            elif basename == "alt_scaffolds":
                # Remember everything in here is alts
                is_alt = True
            else:
                # Otherwise the FASTAs we find will represent primary sequences
                is_alt = False
                
            for item in connection.list_children(scaffold_category):
                # For every tiem in there (FASTA directory, chr2acc,
                # alt_scaffold_placement.txt, etc.)
                
                # Get the base name (last path component)
                item_basename = os.path.basename(item)
                
                Logger.debug("Child item: {}".format(item_basename))
                
                if item_basename == "chr2acc":
                    # This file contains primary chromosome names
                    Logger.info("Chromosome names: {}".format(item))
                    plan.add_primary_scaffold_names(connection.get_url(item))
                elif item_basename == "alt_scaffold_placement.txt":
                    # This file contains alt scaffold placements
                    Logger.info("Alt placements: {}".format(item))
                    plan.add_alt_scaffold_placements(connection.get_url(item))
                elif item_basename == "FASTA":
                    # Look for FASTA files in there
                    for file_path in connection.list_children(item):
                        # For each file
                        if (file_path.endswith(".fa") or
                            file_path.endswith(".fa.gz") or
                            file_path.endswith(".fna") or
                            file_path.endswith(".fna.gz")):
                            
                            # This is a legit-looking FASTA (and probably not a
                            # directory)
                            
                            if is_alt:
                                # Add this FASTA as containing alt scaffolds
                                Logger.info("Alt FASTA: {}".format(file_path))
                                plan.add_alt_scaffold_fasta(connection.get_url(
                                    file_path))
                            else:
                                # Add this FASTA as containing primary scaffolds
                                Logger.info("Primary FASTA: {}".format(
                                    file_path))
                                plan.add_primary_scaffold_fasta(
                                    connection.get_url(file_path))
                                    
    
    # Now we have added everything to the plan that comes from the GRC-format
    # assembly. We still may need VCFs, but that's up to the caller
    return plan
    
def main_job(job, options, plan):
    """
    Root Toil job. Right now does nothing.
    """
    # TODO: implement
    pass
   
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    
    # Build the plan on the head node
    plan = create_plan(options.assembly_structure)
    
    return 1
    
    # Make a root job
    root_job = Job.wrapJobFn(main_job, options, plan,
        cores=1, memory="1G", disk="1G")
    
    # Run it and see how many jobs fail. Automatically handles RealtimeLogger
    # messages
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    return 0
    
def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        
        
        
        
        
        
        
        
        
        

