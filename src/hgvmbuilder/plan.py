#hgvm-builder plan.py: Represent reference graph build plans as objects

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
        
        
    
    def add_primary_scaffold_fasta(self, url):
        """
        Adds the (possibly gzipped) FASTA at the given URL as a FASTA known to
        contain only primary scaffolds (the assembled chromosomes, things like
        chrM, and other top-level unplaced scaffolds).
        
        FASTA records should be named as accession.version.
        """
        
        raise NotImplementedError
        
    def add_primary_scaffold_names(self, url):
        """
        Adds the TSV at the given URL, in GRC chr2acc format (chromosome
        number/X/Y/MT and accession.version, with #-comments) to the database of
        chromosome to accession.version mappings. This database is then used to
        interpret VCFs.
        """
        
        raise NotImplementedError
        
    def add_alt_scaffold_fasta(self, url):
        """
        Adds the (possibly gzipped) FASTA at the given URL as a FASTA known to
        contain only alt scaffolds (alt locus scaffolds and patches).
        
        FASTA records should be named as accession.version.
        
        Alt scaffolds that don't eventually have parents in the files added with
        add_alt_scaffold_placement will cause an error when the plan is
        executed.
        """
        
        raise NotImplementedError
        
    def add_alt_scaffold_placement(self, url):
        """
        Adds the TSV at the given URL, in GRC alt_scaffold_placement.txt format,
        as a database of child alt scaffold to parent primary scaffold mappings.
        Child alt scaffolds will be re-mapped, but this file will restrict re-
        mapping to the appropriate chromosome, so it can be done in parallel.
        """
        
        raise NotImplementedError
        
    def add_variants(self, chromosome_name, url): 
        """
        Adds the gzipped VCF at the given URL as variants to be applied to the
        given chromosome number or X/Y/MT name. Only one VCF is allowed per
        chromosome.
        """
        
        raise NotImplementedError
        
    def add_sample(self, fastq_url):
        """
        Adds the given high coverage sample FASTQ as a sample to be aligned to
        and used to augment the graph.
        
        TODO: Will be implemented after VCFs and assembled contigs work
        """
        
        raise NotImplementedError














    
    
