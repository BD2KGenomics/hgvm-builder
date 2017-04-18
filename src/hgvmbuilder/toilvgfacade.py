# hgvm-builder toilvgfacade.py: Provide a function-argument-based toil-vg API

"""
toil-vg curtrently has lots of cases where low-level functions depend on
command-line arguments in the options object. To make toil-vg targets callable
on arbitrary Toil file IDs, we need wrappers.

To use this facade, run add_options() on your argparse parser before parsing
arguments, initialize() on your resulting options namespace on the master, and
the various _job functions as Toil jobs to actually do stuff.

"""

import os
import os.path
import logging
import urlparse
import shutil
import argparse

import toil_vg.vg_common
import toil_vg.vg_config
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
        
        if name.strip("-") not in self.blacklist:
            # Add it!
            return self.real_parser.add_argument(name, *args, **kwargs)
 
# What options don't we want to pass through to/from the command line? Don't add
# the leading dashes. Holds a dict from toil vg operation type to the options
# that should be removed.
option_blacklist = {
    "wrapperscript": {"out_store", "tool"},
    "common": {"force_outstore"},
    "index": {"index_name", "graphs", "chroms"}
}

# TODO: Un-blacklist --config and add logic to import the config file and send
# it via the file store to the nodes that actually use toil-vg. Or otherwise
# require special prep code to be run on the master to use this library.

def add_options(parser):
    """
    Given an argparse parser or option group, add all the toil-vg configuration
    options (for extra vg command flags, Docker containers, and so on).
    
    """
    
    # Add all the non-blacklisted toil-vg common options
    common_group = parser.add_argument_group("Toil VG configuration",
        "options to configure the Toil VG wrapper")
    toil_vg.vg_common.add_container_tool_parse_args(OptionFilter(common_group,
        option_blacklist["common"]))
    toil_vg.vg_common.add_common_vg_parse_args(OptionFilter(common_group,
        option_blacklist["common"]))
    
    # Add all the non-blacklisted vg index options to this group
    index_group = parser.add_argument_group("VG Indexing",
        "options to configure involations of `vg index`")
    toil_vg.vg_index.index_parse_args(OptionFilter(index_group,
        option_blacklist["index"]))
        
def initialize(options):
    """
    Start up the Toil VG system on the master. Imports a bunch of config file
    defaults into the options.
    
    """
    
    logging.info("Using Toil VG from {}".format(toil_vg.__file__))
    
    # Apply the config file
    processed_options = toil_vg.vg_config.apply_config_file_args(options)
    # Apply the changes back to the original options
    options.__dict__ = processed_options.__dict__
    
    # Make a command runner that uses Docker (or Singularity)
    options.drunner = toil_vg.vg_common.ContainerRunner(
        container_tool_map = toil_vg.vg_common.get_container_tool_map(options))
        
def sanitize_options(cli_options):
    """
    Since Toil VG uses the command line options namespace thingy as a sort of
    general context, we will need to feed one into every call.
    
    However, since we removed some command-line options, our caller might feed
    us an options object with those set (because it re-used those option names).
    So we have to strip them out.
    """
    
    # We'll fill this in
    sanitized = argparse.Namespace()
    
    # We compute a global blacklist of options that some toil vg function
    # shouldn't get.
    global_blacklist = set()
    for local_blacklist in option_blacklist.itervalues():
        for item in local_blacklist:
            # We should strip this out
            global_blacklist.add(item)
    
    for key, value in vars(cli_options).iteritems():
        # For everything we got fed
        if key.strip("-") in global_blacklist:
            # Blacklisted options get skipped
            continue
        
        # Copy everything else over
        vars(sanitized)[key] = value
        
    return sanitized
    
def xg_index_job(job, options, vg_ids):
    """
    Index the given VG graphs into an XG file. Returns the ID of the XG file.
    Automatically sets the correct resource requirements based on the config
    passed via options.
    
    Internally uses toil_vg to perform the indexing.
    """
    
    # Do any options manipulation we need to do
    
    # Strip out stuff we don't want and apply config defaults
    options = sanitize_options(options)
    
    # Add the outstore, which we have sort of disabled. It insists on writing
    # stuff, so just drop it in the current directory. It doesn't read it back.
    options.out_store = "file:."

    # Don't use it instead of the filestore
    options.force_outstore = False
    
    # Pretend we're the pipeline tool
    options.tool = "pipeline"
    
    # Add stuff that toil vg index uses
    
    # options.chroms has to have a name for every graph, to save it under in the
    # local temp dir.
    options.chroms = ["graph{}".format(i) for i in xrange(len(vg_ids))]
    
    # options.index_name has to have the basename for the .xg in the local temp
    # dir.
    options.index_name = "xgindex"
    
    return job.addChildJobFn(run_xg_indexing_wrapper, options,
        vg_ids, cores=options.xg_index_cores, memory=options.xg_index_mem,
        disk=options.xg_index_disk).rv()
        
# Toil explodes when we try to pass it functions from other modules, so just
# wrap all of them for now. See
# <https://github.com/BD2KGenomics/toil/issues/1633>
        
def run_xg_indexing_wrapper(*args, **kwargs):
    return toil_vg.vg_index.run_xg_indexing(*args, **kwargs)
    




















