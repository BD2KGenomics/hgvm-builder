#!/usr/bin/env python
# hgvm-builder chromMatcher.py: match chromosomes from one assembly to another
"""

Takes two assembly FASTAs, a reference and a query, and matches each contig in
the query to one or more contigs in the reference by aligning parts of it.

The reference should have already been indexed with `bwa index`.

"""

import argparse
import logging
import os
import os.path
import re
import shutil
import shlex
import subprocess
import sys
import collections
import tempfile
import StringIO

import Bio.SeqIO
import Bio.Sequencing.Applications

import pysam

import tsv

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
    
    parser.add_argument("ref_name",
        help="name of BWA-indexed reference FASTA")
    parser.add_argument("query_name",
        help="name of query FASTA")
        
    parser.add_argument("--chunk_size", type=int, default=10000,
        help="size of chunks to align")
    parser.add_argument("--match_threshold", type=float, default=0.95,
        help="min score for a hit as a fraction of possible score")
    parser.add_argument("--out_file", type=argparse.FileType("w"),
        default=sys.stdout,
        help="TSV file of contig pairings to write")
    parser.add_argument("--debug", action="store_true",
        help="add extra debugging comments to output TSV")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    writer = tsv.TsvWriter(options.out_file)
    
    for record in Bio.SeqIO.parse(options.query_name, "fasta"):
        # For each record to place
        
        # Where do we put it
        hit_contigs = set()
        
        for chunk_start in xrange(0, len(record.seq), options.chunk_size):
            # Pull out each chunk of the appropriate max size
            chunk = record[chunk_start : chunk_start + options.chunk_size]
            
            # Make a temp file for it
            temp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)
            
            # Save to the temp file
            Bio.SeqIO.write(chunk, temp_file, "fasta")
            temp_file.close()
            
            # Make the BWA command line
            bwasw_cmd = Bio.Sequencing.Applications.BwaBwaswCommandline(
                reference=options.ref_name, read_file=temp_file.name)
            
            # Run it
            child = subprocess.Popen(shlex.split(str(bwasw_cmd)), stdout=subprocess.PIPE)
            
            # Read the SAM file
            alignment_file = pysam.AlignmentFile(child.stdout)
            
            for read in alignment_file.fetch(until_eof=True):
                # For each read, get its contig and score
                score = dict(read.get_tags())["AS"]
                contig = read.reference_name
                
                if options.debug:
                    # Report details on the match
                    writer.comment("{}@{} -> {}: {}".format(record.id,
                        chunk_start, contig, score))
                
                if score > len(read.query_sequence) * options.match_threshold:
                    # Call it a good enough match
                    hit_contigs.add(contig)
                
            alignment_file.close()
            child.wait()
            
            os.unlink(temp_file.name)
            
            if len(hit_contigs) > 0:
                # Only go until the first good hit. Tehn move on to the next
                # record
                break
                
        for ref_contig in hit_contigs:
            # Say each mapping of query contig to ref contig that we found
            writer.line(record.id, ref_contig)
    
    return 0
    
def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        
        
        
        
        
        
        
        
        
        

