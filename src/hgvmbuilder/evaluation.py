#hgvm-builder evaluation.py: Represent reference graph evaluation regimines as objects

import collections
import logging

from .directory import Directory

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
        
        # Then these will hold the IDs
        self.fastq_ids = []
        self.control_graph_id = None
        # This maps from vg ID to xg index ID
        self.xg_ids = {}
        # This maps from vg ID to gcsa and lcp index ID pair
        self.gcsa_lcp_ids = {}
        
        
    def add_fastq(self, fastq_url):
        """
        Add the given FASTQ either as a single FASTQ or as a memebr of the PASTQ
        pair.
        """
        
        assert(len(self.fastq_urls) <= 2)
        self.fastq_urls.append(fastq_url)
        
    def set_eval_sequences(self, url):
        """
        Set the URL from which to obtain the sequences to be realigned to the
        sample graph for evaluating it.
        """
        
        self.eval_sequences_url = url
        
    def set_control_graph(self, url):
        """
        Set the control graph to be compared to the experimental graph.
        """
        
        self.control_graph_url = url
        
    def set_control_graph_xg(self, url):
        """
        Set the XG index for the control graph.
        """
        
        self.control_graph_xg_url = url
        
    def set_control_graph_gcsa(self, gcsa_url, lcp_url):
        """
        Set the GCSA/LCP index for the control graph.
        """
        
        self.control_graph_gcsa_url = gcsa_url
        self.control_graph_gcsa_lcp_url = lcp_url
        
    def get_fastq_ids(self):
        """
        Get the 0-2 FASTQ file IDs to use as a list.
        """
        
        return self.fastq_ids
        
    def get_control_graph_id(self):
        """
        Return the file ID of the control graph, or None if no such graph was
        specified.
        
        """
        
        return self.control_graph_id
        
    def get_xg_id(self, vg_id):
        """
        Return the file ID for the given vg graph ID, or None if no index
        exists.
        """
        
        return self.xg_ids.get(vg_id, None)
        
    def get_gcsa_id(self, vg_id):
        """
        Return the GCSA index ID for the given vg graph ID, or None if no index
        exists.
        """
        
        if self.gcsa_lcp_ids.has_key(vg_id):
            # We have a GCSA/LCP, so return the first ID (GCSA)
            return self.gcsa_lcp_ids[vg_id][0]
        else:
            return None
            
    def get_lcp_id(self, vg_id):
        """
        Return the GCSA LCP index ID for the given vg graph ID, or None if no
        index exists.
        
        """
        
        if self.gcsa_lcp_ids.has_key(vg_id):
            # We have a GCSA/LCP, so return the second ID (LCP)
            return self.gcsa_lcp_ids[vg_id][1]
        else:
            return None
            
    def get_eval_sequences_id(self):
        """
        Return the ID for the evaluation sequences file.
        """
        
        return self.eval_sequences_id
        
    def get_packaged_control(self):
        """
        Package up the control graph from the evaluation plan as a Directory.
        Assumes the control index was already provided.
        
        """
        
        vg_id = self.get_control_graph_id()
        
        package = Directory({
            "hgvm.vg": vg_id,
            "hgvm.xg": self.get_xg_id(vg_id),
            "hgvm.gcsa": self.get_gcsa_id(vg_id),
            "hgvm.gcsa.lcp": self.get_lcp_id(vg_id)
        })
        
        # TODO: Manifest!
        
        return package
        
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
            
        if self.control_graph_url is not None:
            # Grab the control graph since it exists
            self.control_graph_id = import_function(self.control_graph_url)
            
            Logger.info("Imported control graph {} as {}".format(
                self.control_graph_url, self.control_graph_id))
            
            if self.control_graph_xg_url is not None:
                # Grab the control graph xg index since it exists
                self.xg_ids[self.control_graph_id] = import_function(
                    self.control_graph_xg_url)
                
            if (self.control_graph_gcsa_url is not None and
                self.control_graph_gcsa_lcp_url is not None):
                # Grab the control graph gcsa/lcp index since it exists
                self.gcsa_lcp_ids[self.control_graph_id] = (
                    import_function(self.control_graph_gcsa_url),
                    import_function(self.control_graph_gcsa_lcp_url))
        else:
            # No control graph
            self.control_graph_id = None
            
class SVRecallPlan(object):
    """
    Represents a plan to evaluate a graph reference for structural variant
    calling, using a truth-set-based evaluation method. The graph reference is
    used as the basis for variant calling, and then the portion of structural
    variants that the sample is marked as having in the truth set that are
    called as present in the sample is calculated.
    
    Consists of a pair of FASTQ files for the variant calling, and a sample name
    to look up in the original input VCFs to search for called structural
    variants.
    
    Instead of the FASTQs, you can also directly specify an aligned GAM.
    
    """
    
    def __init__(self):
        """
        Make a new empty evaluation plan
        """
        
        # We have the samples we're supposed to be using
        self.sample_names = []
        
        # We have one set of fields that are populated as the plan is built, and
        # another set of fields that get populated when we import the input
        # files into Toil.
        
        # These all hold the URLs that we get given when building up the plan.
        # This will hold lists of up to 2 FASTQs by sample name
        self.fastq_urls = collections.defaultdict(list)
        # This will hold single GAM URLs by sample name
        self.gam_urls = {}
        # This will hold single call VCF URLs by sample name
        self.vcf_urls ={}
        
        # Then these will hold the IDs, following the same shape.
        # For the FASTQs
        self.fastq_ids = collections.defaultdict(list)
        # And the GAMs
        self.gam_ids = {}
        # And the VCFs
        self.vcf_ids = {}
                    
        
    def add_sample_fastq(self, sample_name, fastq_url):
        """
        Add the given FASTQ either as a single FASTQ or as a memebr of the FASTQ
        pair for a sample.
        """
        
        assert(len(self.fastq_urls[sample_name]) < 2)
        self.fastq_urls[sample_name].append(fastq_url)
        
    def add_sample_gam(self, sample_name, gam_url):
        """
        Set the given GAM file to use instead of the FASTQs for the given
        sample.
        
        """
        
        self.gam_urls[sample_name] = gam_url
        
    def add_sample_vcf(self, sample_name, vcf_url):
        """
        Set the given VCF file to use instead of the FASTQsor GAM for the given
        sample.
        
        """
        
        self.vcf_urls[sample_name] = vcf_url
        
    def add_sample_name(self, name):
        """
        Add a sample that we will call and evaluate SVs on.
        """
        
        self.sample_names.append(name)
        
    def get_fastq_ids(self, sample_name):
        """
        Get the 0-2 FASTQ file IDs to use for this sample as a list.
        """
        
        return self.fastq_ids[sample_name]
        
    def get_gam_id(self, sample_name):
        """
        Get the GAM file ID to use for the given sample, or None.
        """
        
        return self.gam_ids.get(sample_name, None)
        
    def get_vcf_id(self, sample_name):
        """
        Get the VCF file ID to use for the given sample, or None.
        """
        
        return self.vcf_ids.get(sample_name, None)
        
    def get_sample_names(self):
        """
        Return the sample names we are using to evaluate SV calling.
        
        """
        
        return self.sample_names
        
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
        
        for sample, urls in self.fastq_urls.iteritems(): 
            # Grab the FASTQs
            self.fastq_ids[sample] = [import_function(url) for url in urls]
        
        for sample, url in self.gam_urls.iteritems():
            # Grab the GAM
            self.gam_ids[sample] = import_function(url)
            
        for sample, url in self.vcf_urls.iteritems():
            # Grab the VCF
            self.vcf_ids[sample] = import_function(url)
        
    
