#hgvm-builder grcparser.py: utilities for parsing GRC-format assembly structure

import logging
import urlparse

from .ftputil import FTPOrFilesystemConnection

# Get a submodule-global logger
Logger = logging.getLogger("grcparser")

def parse(plan, assembly_structure):
    """
    Given a plan.ReferencePlan to fill in and a GRC assembly structure format
    URL like ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.24_GRCh38.p9/G
    CA_000001405.24_GRCh38.p9_assembly_structure to parse, connect to the URL
    and fill in the plan with info from the assembly.
    """
    
    # Open the URL
    connection = FTPOrFilesystemConnection(assembly_structure)
    
    for assembly_unit_root in connection.list_children(""):
        # For each assembly unit
        Logger.debug("Assembly unit: {}".format(assembly_unit_root))
        for scaffold_category in connection.list_children(assembly_unit_root):
            # For everything under that
            
            # Get the base name (last path component)
            basename = scaffold_category.split("/")[-1]
            
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
                item_basename = item.split("/")[-1]
                
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
