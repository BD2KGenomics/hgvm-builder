# hgvm-builder toilvgfacade.py: Provide a function-argument-based toil-vg API

# toil-vg curtrently has lots of cases where low-level functions depend on
# command-line arguments in the options object. To make toil-vg targets callable
# on arbitrary Toil file IDs, we need wrappers.

import os
import os.path
import logging
import urlparse
import shutil

import toil_vg.vg_common
import toil_vg.vg_index


Logger = logging.getLogger("toilvgfacade")

class OptionFilter(object):
    """
    Can wrap an ArgumentParser or other such class and drop options on a
    blacklist/accept only options on a whitelist.
    
    """
    
    def __init__(self, real_parser, blacklist=[]):
        """
        Wrap the given actual parser with an add_option method.
        """
        
        # Save the parser
        self.real_parser = real_parser
        
        # Save the blacklist
        self.blacklist = set(blacklist)
        
    def add_argument(self, name, *args, **kwargs):
        """
        Add the given argument, if its name passes the filters.
        
        """
        
        if name not in self.blacklist:
            # Add it!
            return self.real_parser.add_argument(name, *args, **kwargs)

def add_options(parser):
    """
    Given an argparse parser or option group, add all the toil-vg configuration
    options (for extra vg command flags, Docker containers, and so on).
    
    """
    
    # Add all the toil-vg common options
    common_group = parser.add_argument_group("Toil VG configuration",
        "options to configure the Toil VG wrapper")
    toil_vg.vg_common.add_container_tool_parse_args(common_group)
    toil_vg.vg_common.add_common_vg_parse_args(common_group)
    
    # Add all the vg index options to this group
    index_group = parser.add_argument_group("VG Indexing",
        "options to configure involations of `vg index`")
    # Except for these
    index_blacklist = {"--index_name", "--graphs", "--chroms"}
    toil_vg.vg_index.index_parse_args(OptionFilter(index_group,
        index_blacklist))




















