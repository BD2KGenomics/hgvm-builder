#!/usr/bin/env python
# hgvm-builder importsam.py: Command-line tool to download SAM/BAM/CRAM reads
"""

Takes one or more SAM/BAM/CRAM URLs and contig names, and extracts the reads on
those contigs and with pair partners on those contigs.

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
import shutil
import urlparse

import tsv

import toil
import toil.version
from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from .directory import Directory
from .toilpromise import ToilPromise
from . import toilvgfacade

from . import smartSam2Fastq

# Get a submodule-global logger
Logger = logging.getLogger("importsam")

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
    
    parser.add_argument("--contig", default=None,
        help="download just pairs touching this contig")
    
    parser.add_argument("--sam_url", default=[], action="append",
        help="URL to download reads from, accessible from Toil nodes")
        
    # Output
    parser.add_argument("out_url",
        help="file: or other Toil-supported URL to place results in")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
def extract_job(job, options, sam_url):
    """
    Extract and fix up the given SAM/BAM/CRAM reads by URL.
    
    Return a pair of FASTQ file IDs.
    """
    
    # Let's just download the whole bam
    sorted_bam = job.fileStore.getLocalTempDir() + "/sorted.bam"
    
    # We need a prefix for temp files
    temp_prefix = sorted_bam + ".part"
    
    RealtimeLogger.info("Sort {} to BAM".format(sam_url))
    
    # Sort reads by name to a BAM file. If we don't give a temp file prefix it
    # tries to write the temp files back to the FTP.
    options.drunner.call(job, [["samtools", "sort", "-n", "-o",
        sorted_bam, "-T", temp_prefix, sam_url]])
        
    # Save to file store
    bam_id = job.fileStore.writeGlobalFile(sorted_bam)
    
    # Convert and return FASTQs
    return job.addChildJobFn(convert_job, options, sam_url, bam_id,
        cores=4, memory="16G", disk="1000G")
        
def convert_job(job, options, sam_url, bam_id):
    """
    Subset and convert BAM to FASTQ pair. Returns FASTQ IDs.
    """
    
    # Read the BAM back
    sorted_bam = job.fileStore.readGlobalFile(bam_id)
        
    RealtimeLogger.info("Subset {} to SAM".format(sam_url))
        
    # Then stream to SAM and select just the reads we want
    subset_sam = job.fileStore.getLocalTempDir() + "/subset.sam"
    
    # We start out with just a view pipeline
    sam_command = [["samtools", "view", sorted_bam]]
    
    if options.contig is not None:
        # Subset to this contig and related alts with awk
        sam_command.append(["awk",
            ("{if ($3 ~ /" + options.contig + "(_.*)?$/ || $7 ~ /" + 
            options.contig + "(_.*)?$/) print}")])
    
    options.drunner.call(job, sam_command, outfile=open(subset_sam, "w"))
        
    RealtimeLogger.info("Convert {} to FASTQ".format(sam_url))
    
    with job.fileStore.writeGlobalFileStream() as (fq1_handle, fq1_id):
        with job.fileStore.writeGlobalFileStream() as (fq2_handle, fq2_id):
    
            # Then prep options for running the converter script in this Python
            convert_options = argparse.Namespace()
            convert_options.input_sam = open(subset_sam, "r")
            convert_options.fq1 = fq1_handle
            convert_options.fq2 = fq2_handle
            convert_options.drop_secondary = True
            convert_options.expect_paired = True
            convert_options.interleaved = False
        
            smartSam2Fastq.run(convert_options)
            
            return fq1_id, fq2_id
    
   
def main_job(job, options, sam_urls):
    """
    Root job of the Toil workflow. Download the sample URLs.
    
    Returns a Directory containing a bunch of output files.
    
    """
    
    RealtimeLogger.info("Main job starting")
    
    # Make sure we can use samtools
    options.drunner.call(job, [["samtools", "--version"]])
    
    # We'll fill this with promises for subdirectories by sample filename
    subdir_promises = {}
    
    for sam_url in sam_urls:
        # Work out the base filename
        sam_filename = os.path.basename(urlparse.urlparse(sam_url).path)
        
        # Go download and convert the reads, and stick the FASTQs in a directory
        subdir_promises[sam_filename] = ToilPromise.wrap(
            job.addChildJobFn(extract_job, options, sam_url,
            cores=1, memory="16G", disk="1000G")
        ).then(lambda (fq1, fq2): Directory({
            "fq1.fastq": fq1, "fq2.fastq": fq2}))
        
    # Mount each subdirectory under its original sam/bam/cram filename
    return ToilPromise.all(subdir_promises
        ).then(lambda dirs: Directory().mount_all(dirs)
        ).unwrap_result()
        
    
    
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
        
            # Don't import on the master. Let the nodes handle the download.
        
            # Make a root job
            root_job = Job.wrapJobFn(main_job, options, options.sam_url,
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
        
        
        
        
        
        
        
        
        
        

