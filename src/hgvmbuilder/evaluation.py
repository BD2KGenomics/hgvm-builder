#hgvm-builder evaluation.py: Represent reference graph evaluation regimines as objects

import logging
import urllib2
import collections

import tsv

Logger = logging.getLogger("evaluation")

class EvaluationPlan(object):
    """
    Represents a plan to evaluate a graph reference, using a truth-set-free
    evaluation method. The graph reference is used as the basis for variant
    calling, and then the sample graph produced from the calling process is used
    as an alignment target for other sequences from the same sample. The graph
    is evaluated based on how well it embeds the test sequences.
    
    Consists of a pair of FASTQ files for the variant calling, a file of
    evaluation sequences for realigning, and a control graph with indexes that
    the graph reference is compared against.
    
    """
    
    def __init__(self):
        """
        Make a new empty evaluation plan
        """
        
        # We have one set of fields that are populated as the plan is built, and
        # another set of fields that get populated when we import the input
        # files into Toil.
        
        # These all hold the URLs that we get given when building up the plan
        # This will hold 2 FASTQ URLs.
        self.fastq_urls = []
        
        # This holds the URL to the evaluation sequence file
        self.eval_sequences_url = None
        
        # This holds the control graph URL to the vg file
        self.control_graph_url = None
        # This holds the URL to its XG file
        self.control_graph_xg_url = None
        # And to its GCSA file
        self.control_graph_gcsa_url = None
        # And to its GCSA LCP file
        self.control_graph_gcsa_lcp_url = None
        
    def bake(self, import_function):
        """
        "Bake" the plan by importing data files into a file storage system.
        Requires a function that can take a URL and load it into a file storage
        system, returning an ID from which the file can be retrieved. Generally
        you would get this from something like:
        
            with toil.common.Toil(options) as toil_instance:
                plan.bake(lambda url: toil_instance.importFile(url))
        
        But we don't want to attach directly to Toil here, so you can pass in
        anything you want.
        """
        
        # Import everything
        
        # Grab the FASTQs
        self.fastq_ids = [import_function(url) for url in self.fastq_urls]
        
        if self.eval_sequences_url is not None:
            # Grab the eval sequences since they exist
            self.eval_sequences_id = import_function(self.eval_sequences_url)
        else:
            self.eval_sequences_id = None
            
        # TODO: finish this!
                   
        
        
        














    
    
