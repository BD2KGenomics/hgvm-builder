#hgvm-builder thousandgenomesparser.py: utilities for parsing VCF directories

import logging
import re

from .ftputil import FTPOrFilesystemConnection

# Get a submodule-global logger
Logger = logging.getLogger("thousandgenomesparser")

def parse(plan, vcf_root):
    """
    Given a plan.ReferencePlan to fill in and a URL to a directory of VCFs like 
    http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38
    _positions, list all the VCFs and add them to the plan.
    
    Assumes that each VCF has "chr[0-9A-Za-z]+" somewhere in its name,
    identifying the chromosome it belongs to.
    
    """
    
    # Open the URL
    connection = FTPOrFilesystemConnection(vcf_root)
    
    Logger.info("Connected to {}".format(vcf_root))
    
    # Define the chromosome identification regex
    chrom_finder = re.compile("chr([0-9A-Za-z]+)")
    
    for item in connection.list_children(""):
        # For each file that might be a VCF
        
        if (not item.endswith(".vcf")) and (not item.endswith(".vcf.gz")):
            # Not a VCF
            continue
            
        # Find a match for the chromosome name pattern, if any exists            
        match = chrom_finder.search(item)
        if match is not None:
            # We found one. Pull out the chromosome name without chr
            name = match.group(1)
            
            Logger.info("Chromosome {} VCF: {}".format(name, item))
            
            # Add the VCF to the plan
            plan.add_variants(name, connection.get_url(item))
