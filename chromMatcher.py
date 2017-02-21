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
    parser.add_argument("--max_children", type=int, default=10,
        help="number of bwa children to run")
    parser.add_argument("--batch_size", type=int, default=1000,
        help="number of chunks to align at once")
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
    
    
def run_batch(options, batch_state, matches, writer):
    """
    Takes a list of (SeqRecord, offset) pairs and a defaultdict form sequence ID
    to set of hits.
    
    Run a batch with the latest piece form every state. Throw out the states
    that finish. Returns the updated list of (SeqRecord, new offset) pairs for
    records that still have no hits.
    
    Writes output to the given TsvWriter.
    """
    
    # Break up batches among children
    child_batches = [[] for i in xrange(options.max_children)]
    for i, batch in enumerate(batch_state):
        child_batches[i % options.max_children].append(batch)
        
    # This holds the child subprocesses
    children = []
    
    # This holds the temp files to clean up
    temp_files = []
    
    for child_batch_state in child_batches:
    
        # Make a temp file for the FASTA to align
        temp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)
        temp_files.append(temp_file)
        
        for (record, offset) in child_batch_state:
            # Pull out each chunk of the appropriate max size
            chunk = record[offset : offset + options.chunk_size]
            
            # Save to the temp file
            Bio.SeqIO.write(chunk, temp_file, "fasta")
            
        temp_file.close()
        
        # Make the BWA command line
        bwa_cmd = Bio.Sequencing.Applications.BwaBwaswCommandline(
            reference=options.ref_name, read_file=temp_file.name)
        
        # Run it
        child = subprocess.Popen(shlex.split(str(bwa_cmd)), stdout=subprocess.PIPE)
        children.append(child)
        
    for child in children:
        # Now handle all the children one at a time  
        
        # Read the SAM file
        alignment_file = pysam.AlignmentFile(child.stdout)
        
        for read in alignment_file.fetch(until_eof=True):
            # For each read, get its contig and score
            score = dict(read.get_tags())["AS"]
            contig = read.reference_name
            
            threshold_score = len(read.query_sequence) * options.match_threshold
            
            if options.debug:
                # Report details on the match
                writer.comment("{} -> {}: {}/{}".format(read.query_name, contig,
                    score, threshold_score))
            
            if score >= threshold_score:
                # Call it a good enough match
                matches[read.query_name].add(contig)
        
        alignment_file.close()
        child.wait()
        
    for temp_file in temp_files:
        # Clean up temp files
        os.unlink(temp_file.name)
    
    # Now update the state
    new_batch_state = []
    for (record, offset) in batch_state:
        if len(matches[record.id]) > 0:
            # We have a hit for this record, so we don't need to run it any
            # more
            for ref_contig in matches[record.id]:
                # Say each mapping of query contig to ref contig that we
                # found
                writer.line(record.id, ref_contig)
        else:
            # No hits, so we need to keep it
            # Where will we align next?
            new_offset = offset + options.chunk_size
            
            if new_offset < len(record.seq):
                new_batch_state.append((record, new_offset))
                
            else:
                # No hits before end of contig
                writer.comment("Record {} has no hits".format(
                    record.id))
    return new_batch_state
   
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    writer = tsv.TsvWriter(options.out_file)
    
    # We use this to keep our state: pairs of SeqRecord and current offset
    batch_state = []
    
    # We use this to keep track of the matches found for each query sequence ID
    matches = collections.defaultdict(set)
    
    # OK so we need to run that update function a bunch and feed in data
    
    for record in Bio.SeqIO.parse(options.query_name, "fasta"):
        # For each record to place
        
        while len(batch_state) >= options.batch_size:
            # Run through existing batches until there's room
            batch_state = run_batch(options, batch_state, matches, writer)
            
        # Start a new fake thread on this record
        batch_state.append((record, 0))
        
    # Now run all the records to completion
    while len(batch_state) > 0:
        batch_state = run_batch(options, batch_state, matches, writer)
        
    return 0
    
def entrypoint():
    """
    0-argument entry point for setuptools to call.
    """
    
    # Provide main with its arguments and handle exit codes
    sys.exit(main(sys.argv))
    
if __name__ == "__main__" :
    entrypoint()
        
        
        
        
        
        
        
        
        
        

