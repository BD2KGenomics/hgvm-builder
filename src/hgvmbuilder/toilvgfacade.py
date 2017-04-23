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
import timeit

import toil_vg.vg_common
import toil_vg.vg_config
import toil_vg.vg_index

from toil.realtimeLogger import RealtimeLogger

from .toilpromise import ToilPromise


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

    # Don't use outstore instead of the file store
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
    
    return job.addChildJobFn(toil_vg.vg_index.run_xg_indexing, options,
        vg_ids, cores=options.xg_index_cores, memory=options.xg_index_mem,
        disk=options.xg_index_disk).rv()
        
def gcsa_index_job(job, options, vg_ids, primary_path_names=None):
    """
    Index the given graphs into a GCSA/LCP index, and return a pair of file IDs
    for the GCSA and the LCP files.
    
    Will prune the graph before indexing unless options.prune_opts is explicitly
    set as an empty list.
    
    """
    
    # Do any options manipulation we need to do
    
    # Strip out stuff we don't want and apply config defaults
    options = sanitize_options(options)
    
    # Add the outstore, which we have sort of disabled. It insists on writing
    # stuff, so just drop it in the current directory. It doesn't read it back.
    options.out_store = "file:."

    # Don't use outstore instead of the file store
    options.force_outstore = False
    
    # Pretend we're the pipeline tool
    options.tool = "pipeline"
    
    # Add stuff that toil vg index uses
    
    # options.graphs has to have a name for every graph, to save it under in the
    # local temp dir.
    options.graphs = ["graph{}".format(i) for i in xrange(len(vg_ids))]
    
    # We also need a "chroms" giving the primary path for each graph. It's OK if
    # the path doesn't exist in a given graph, but if it does it will be added
    # to the index.
    if primary_path_names is not None:
        # We have primary path names to use
        assert(len(primary_path_names) == len(vg_ids))
        options.chroms = primary_path_names
    else:
        # Fake path names
        options.chroms = ["" for x in vg_ids]
    
    # options.index_name has to have the basename for the .gcsa in the local
    # temp dir.
    options.index_name = "gcsaindex"
    
    return job.addChildJobFn(toil_vg.vg_index.run_gcsa_prep, options, vg_ids,
        cores=options.misc_cores, memory=options.misc_mem,
        disk=options.misc_disk).rv()
        
def id_range_job(job, options, vg_id):
    """
    Find the first and last ID in the given VG file and return them as a tuple.
    
    """
    
    # Strip out stuff we don't want and apply config defaults
    options = sanitize_options(options)
    
    # We need an options.chroms, because the job we're running returns an entry
    # form it.
    options.chroms = [None]
    
    # Don't use outstore instead of the file store
    options.force_outstore = False
    
    # Run the job and return the start and end IDs as a pair of ints (dropping
    # the chrom name)
    return ToilPromise.wrap(job.addChildJobFn(toil_vg.vg_index.run_id_range,
        options, 0, vg_id,
        cores=options.misc_cores, memory=options.misc_mem,
        disk=options.misc_disk)
        ).then(lambda (name, start, end): (int(start), int(end))
        ).unwrap_result()
        
def id_increment_job(job, options, vg_id, distance):
    """
    Increment all the node IDs in the given vg graph by the given distance.
    Return a new vg graph file ID.
    
    Not actually in toil-vg, but we put it here so all the vg-touching functions
    can live in one place.
    
    """
    
    if distance == 0:
        # No need to shift at all
        return vg_id
    
    # Strip out stuff we don't want and apply config defaults
    options = sanitize_options(options)
    
    # We need an options.chroms, because the job we're running uses it for local
    # filenames
    options.chroms = ["increment"]
    
    # Don't use outstore instead of the file store
    options.force_outstore = False
    
    return job.addChildJobFn(run_id_increment, options, 0, vg_id, distance,
        cores=options.misc_cores, memory=options.misc_mem,
        disk=options.misc_disk).rv()
    
def run_id_increment(job, options, graph_i, graph_id, distance):
    """
    Actually do the ID incrementing. Is a separate, toil-vg-style job so it
    can be added to toil-vg and so we can set the correct resource requirements.
    
    """
    
    RealtimeLogger.info("Starting graph shift...")
    start_time = timeit.default_timer()
    
    work_dir = job.fileStore.getLocalTempDir()

    # download graph
    graph_filename = os.path.join(work_dir, '{}.vg'.format(
        options.chroms[graph_i]))
    toil_vg.vg_common.read_from_store(job, options, graph_id, graph_filename)

    # Output
    output_graph_filename = graph_filename + '.shifted.vg'
   
    RealtimeLogger.info("Moving {} up by {} to {}".format(
        graph_filename, distance, output_graph_filename))

    with open(output_graph_filename, "w") as out_file:
        command = ['vg', 'ids', '--increment', str(distance),
            os.path.basename(graph_filename)]
        options.drunner.call(job, command, work_dir=work_dir, outfile=out_file)

    # Back to store
    output_graph_id = toil_vg.vg_common.write_to_store(job, options,
        output_graph_filename)
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished graph shift. Process took {} seconds.".format(
        run_time))

    return output_graph_id
    
    
    
    
    
    


















