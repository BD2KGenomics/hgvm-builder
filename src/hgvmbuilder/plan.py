#hgvm-builder plan.py: Represent reference graph build plans as objects

import urllib2

import tsv

class ReferencePlan:
    """
    Represents a plan to build a graph reference, or (Human) Genome Variation
    Map.
    
    Consists of a set of base contigs from an assembly (the primary path for
    each chromosome, plus unplaced things), some variants to apply to those
    contigs, some additional contigs to merge into each of those contigs, and
    some new samples to map in and extract variants from.
    
    Also keeps track of what VCF contigs/chromosomes map to what graph-space
    path names (which should all be in accession.version format), and what alt
    contigs are children of what primary contigs (and so should be aligned to
    them).
    
    After making a plan, and adding FASTAs, VCFs, and metadata files to it with
    the appropriate methods, call bake() to download the metadata and import the
    data files into a Toil filestore.
    
    """
    
    def __init__(self):
        """
        Make a new empty reference build plan
        """
        
        # We have one set of fields that are populated as the plan is built, and
        # another set of fields that get populated when we read the database/ID
        # conversion files and import the FASTA and VCF files into Toil.
        
        # These all hold the URLs that we get given when building up the plan
        
        # This holds possibly-gzipped FASTAs containing primary scaffolds, as a
        # list of URLs
        self.primary_fasta_urls = []
        # This holds a list of the primary scaffold chromosome name/number
        # assignment file URLs
        self.primary_name_urls = []
        # This holds possibly-gzipped FASTAs containing alt scaffolds and
        # patches, as a list of URLs
        self.alt_fasta_urls = []
        # This holds a list of alt scaffold placement database file URLs
        self.alt_placement_urls = []
        # This holds a dict from chromosome name/number to relevant VCF URL.
        self.vcf_urls = {}
        
        # After we load all the databases and import all the input files into
        # Toil, we populate these fields.
        
        # This dict maps from alt scaffold accession.version string to parent
        # scaffold accession.version string
        self.alt_parents = {}
        # This dict maps primary scaffold accession.version to chromosome name,
        # if any.
        self.primary_names = {}
        # This dict maps from chromosome name to a Toil file ID for the VCF for
        # that chromosome, if any.
        self.vcf_ids = {}
        # This list contains the Toil file store IDs for all the input possibly-
        # gzipped FASTAs with primary scaffolds.
        self.primary_ids = []
        # This list contains the Toil file store IDs for all the input possibly-
        # gzipped FASTAs with alt scaffolds.
        self.alt_ids = []
        
        # We're going to have to split the primary FASTAs up and shuffle them to
        # be with the right alt FASTA sequences and VCFs. But we'll do that
        # later, not in the plan setup, because we have to break open a large
        # number of FASTAs.
        
    
    def add_primary_scaffold_fasta(self, url):
        """
        Adds the (possibly gzipped) FASTA at the given URL as a FASTA known to
        contain only primary scaffolds (the assembled chromosomes, things like
        chrM, and other top-level unplaced scaffolds).
        
        FASTA records should be named as accession.version.
        """
        
        self.primary_fasta_urls.append(url)
        
    def add_primary_scaffold_names(self, url):
        """
        Adds the TSV at the given URL, in GRC chr2acc format (chromosome
        number/X/Y/MT and accession.version, with #-comments) to the database of
        chromosome to accession.version mappings. This database is then used to
        interpret VCFs.
        """
        
        self.primary_name_urls.append(url)
        
    def add_alt_scaffold_fasta(self, url):
        """
        Adds the (possibly gzipped) FASTA at the given URL as a FASTA known to
        contain only alt scaffolds (alt locus scaffolds and patches).
        
        FASTA records should be named as accession.version.
        
        Alt scaffolds that don't eventually have parents in the files added with
        add_alt_scaffold_placement will cause an error when the plan is
        executed.
        """
        
        self.alt_fasta_urls.append(url)
        
    def add_alt_scaffold_placement(self, url):
        """
        Adds the TSV at the given URL, in GRC alt_scaffold_placement.txt format,
        as a database of child alt scaffold to parent primary scaffold mappings.
        Child alt scaffolds will be re-mapped, but this file will restrict re-
        mapping to the appropriate chromosome, so it can be done in parallel.
        """
        
        self.alt_placement_urls.append(url)
        
    def add_variants(self, chromosome_name, url): 
        """
        Adds the gzipped VCF at the given URL as variants to be applied to the
        given chromosome number or X/Y/MT name. Only one VCF is allowed per
        chromosome.
        """
        
        if self.vcf_urls.has_key(chromosome_name):
            # We already have a VCF for this chromosome
            raise RuntimeError(
                "Duplicate chromosome name {} adding VCF {} to plan".format(
                chromosome_name, url))
                
        self.vcf_urls[chromosome_name] = url
        
    def add_sample(self, fastq_url):
        """
        Adds the given high coverage sample FASTQ as a sample to be aligned to
        and used to augment the graph.
        
        TODO: Will be implemented after VCFs and assembled contigs work
        """
        
        raise NotImplementedError
        
    def bake(self, import_function):
        """
        "Bake" the plan by downloading metadata and importing data files into a
        file storage system. Requires a function that can take a URL and load it
        into a file storage system, returning an ID from which the file can be
        retrieved. Generally you would get this from something like:
        
            with toil.common.Toil as file_importer:
                plan.bake(lambda url: file_importer.importFile(url))
        
        But we don't want to attach directly to Toil here, so you can pass in
        anything you want.
        """
        
        # Download the metadata files
        
        for name_url in self.primary_name_urls:
            # Download and parse all the primary scaffold chromosome names
            self.parse_name_stream(urllib2.urlopen(name_url))
            
        for placement_url in self.alt_placement_urls:
            # Download and parse all the alt scaffold parent information
            self.parse_placement_stream(urllib2.urlopen(placement_url))
        
        # Import all the data files
        
        for primary_url in self.primary_fasta_urls:
            # Import all the primary FASTAs
            self.primary_ids.append(import_function(primary_url))
            
        for alt_url in self.alt_fasta_urls:
            # Import all the alt FASTAs
            self.alt_ids.append(import_function(alt_url))
        
        
    def set_primary_name(self, accession, name):
        """
        Set the chromosome name for a primary scaffold. Takes the
        accession.version string (like "LOL1234.1") and the chromosome name
        string (like "1" or "MT").
        """

        self.primary_names[accession] = name
        
    def set_alt_parent(self, alt_accession, parent_accession):
        """
        Takes the accession.version string (like "LOL1234.1") fort an alt
        scaffold and a primary scaffold, and makes the primary scaffold the
        parent of the alt scaffold. This means that the alt scaffold will be
        aligned against and merged into the subgraph derived from the primary
        scaffold and ist other alt scaffolds.
        """

        self.alt_parents[alt_accession] = parent_accession
        
    def parse_name_stream(self, stream):
        """
        Parse GRC chr2acc format (TSV of name and accession.version) on the
        given input stream, and make the appropriate primary scaffold name
        assignments.
        """
        
        # Make a TSV reader
        reader = tsv.TsvReader(stream)
        
        for name, accession in reader:
            # Apply each name/accession mapping
            self.set_primary_name(accession, name)
            
    def parse_placement_stream(self, stream):
        """
        Parse GRC alt_scaffold_placement.txt format (TSV with alt
        accession.version on column 4 and parent accession.version on column 7)
        on the given input stream, and make the appropriate parent assignments.
        """
        
        # Make a TSV reader
        reader = tsv.TsvReader(stream)
        
        for parts in reader:
            # Look at each non-comment line
            
            if len(parts) < 7:
                # We can't pull out the parent
                raise RuntimeError(
                    "Insufficient columns in alt scaffold location data")
            
            # Make the alt (1-based column 4) a child of the parent (1-based
            # column 7)
            self.set_alt_parent(parts[3], parts[6])
            
                   
        
        
        














    
    
