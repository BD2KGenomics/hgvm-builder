#hgvm-builder plan.py: Represent reference graph build plans as objects

import logging
import urllib2
import collections

import tsv

Logger = logging.getLogger("plan")

class ReferencePlan(object):
    """
    Represents a plan to build a graph reference, or (Human) Genome Variation
    Map.
    
    Consists of a set of a base graph in VG or HAL format and some variants to
    apply (and the assembly used to interpret those variants).
    
    Also keeps track of what VCF contigs/chromosomes map to what graph-space
    path names (which should all be in accession.version format), and what alt
    contigs are children of what primary contigs (and so should be aligned to
    them).
    
    After making a plan, and adding HALs, VGs, FASTAs, VCFs, and metadata files
    to it with the appropriate methods, call bake() to download the metadata and
    import the data files into a Toil filestore.
    
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
        # This holds a dict from chromosome name/number to a list of relevant
        # VCF URLs.
        self.vcf_urls = collections.defaultdict(list)
        # This holds a dict from VCF URL to VCF index URL
        self.vcf_index_urls = {}
        
        # After we load all the databases and import all the input files into
        # Toil, we populate these fields.
        
        # This dict maps from alt scaffold accession.version string to parent
        # scaffold accession.version string
        self.alt_parents = {}
        # This dict maps primary scaffold accession.version to chromosome name,
        # if any.
        self.accession_to_name = {}
        # This does the reverse
        self.name_to_accession = {}
        # This dict maps from chromosome name to a list of Toil file IDs for the
        # VCFs for that chromosome, if any.
        self.vcf_ids = collections.defaultdict(list)
        # This holds the list of raw, possibly compressed FASTA file IDs for
        # FASTAs with primary sequences.
        self.primary_ids = []
        # This holds the list of raw, possibly compressed FASTA file IDs for
        # FASTAs with alt sequences.
        self.alt_ids = []
        
        # This maps from primary scaffold name to file ID for its uncompressed
        # FASTA. The FASTA may contain other sequences.
        self.uncompressed_primary_fastas = {}
        
        # This maps from alt scaffold name to file ID for its uncompressed
        # FASTA. The FASTA may contain other sequences.
        self.uncompressed_primary_fastas = {}
        
        # This maps from VCF or FASTA file ID to index (.tbi or .fai) file ID
        self.index_ids = {}
        
        # We're going to have to split the primary FASTAs up and shuffle them to
        # be with the right alt FASTA sequences and VCFs. But we'll do that
        # later, not in the plan setup, because we have to break open a large
        # number of FASTAs.
        
        # We also hold files that can bring in whole graphs at various stages of
        # the process
        
        # These are the URLs of HAL files representing graphs to start with
        self.hal_urls = []
        # And these are their IDs after uploading
        self.hal_ids = []
        
        # These are the URLs of VG graphs that we pull in as already converted
        # from HAL
        self.base_vg_urls = []
        # And these are their IDs after uploading
        self.base_vg_ids = []
        
    
    def add_hal(self, url):
        """
        Add a HAL file to base the initial graph on.
        
        """
        
        self.hal_urls.append(url)
        
    def for_each_hal(self):
        """
        In a baked plan, return an iterator over HAL IDs.
        """
        
        return iter(self.hal_ids)
        
    def add_base_vg(self, url):
        """
        Add a VG file to base the initial graph on.
        
        """
        
        self.base_vg_urls.append(url)
        
    def for_each_base_vg(self):
        """
        In a baked plan, return an iterator over base VG IDs.
        """
        
        return iter(self.base_vg_ids)
    
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
        
    def add_alt_scaffold_placements(self, url):
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
        given chromosome number or X/Y/MT name.
        """
        
        self.vcf_urls[chromosome_name].append(url)
        
    def add_variants_index(self, url): 
        """
        Adds the .tbi index at the given URL as the index of variants to be
        applied. Only one index is allowed per VCF, and the index URL must be
        the URL of the VCF with ".tbi" appended. The index may be added before
        or after the VCF itself.
        """
        
        if not url.lower().endswith(".tbi"):
            # Make sure the last 4 characters are what we expect
            raise RuntimeError("Index does not appear to be Tabix: " + url)        
        
        # Drop the last 4 characters to get the VCF name
        self.vcf_index_urls[url[:-4]] = url
        
    def add_sample(self, fastq_url):
        """
        Adds the given high coverage sample FASTQ as a sample to be aligned to
        and used to augment the graph.
        
        TODO: Will be implemented after VCFs and assembled contigs work
        """
        
        raise NotImplementedError
        
    def chromosome_name_to_accession(self, chromosome_name):
        """
        Given a chromosome name, return the corresponding accession for the
        primary contig.
        """
        # Return the changed name, or the unmodified value if there's no mapping
        return self.name_to_accession.get(chromosome_name, chromosome_name)
        
    def accession_to_chromosome_name(self, accession):
        """
        Given an accession, return the chromosome name if any is associated, or
        the accession otherwise.
        """
        # Return the changed name, or the unmodified value if there's no mapping
        return self.accession_to_name.get(accession, accession)

    def for_each_primary_fasta_id(self):
        """
        Iterate over primary input FASTA IDs.
        """
        return iter(self.primary_ids)
        
    def for_each_alt_fasta_id(self):
        """
        Iterate over alt input FASTA IDs.
        """
        return iter(self.alt_ids) 
            
    def for_each_vcf_id_by_chromosome(self):
        """
        Iterate through (chromosome name, VCF file ID) pairs. Each chromosome
        can appear multiple times.
        """
        
        for (chromosome, vcfs) in self.vcf_ids.iteritems():
            for vcf in vcfs:
                # Pair each chromosome with each of its VCF file IDs
                yield (chromosome, vcf)
        
    def for_each_chromosome(self):
        """
        Iterate through chromosome names that have associated VCFs.
        """
        
        return self.vcf_ids.iterkeys()
        
    def get_vcf_ids(self, chrom_name):
        """
        Get the VCF file IDs for the given chromosome name, or [] if no VCFs are
        associated with a name.
        """
        
        return self.vcf_ids[chrom_name]
        
    def get_index_id(self, file_id):
        """
        Given a VCF or FASTA file ID, gets the file ID for the index for that
        file.
        """
        
        return self.index_ids.get(file_id, None)
        
    def bake(self, import_function):
        """
        "Bake" the plan by downloading metadata and importing data files into a
        file storage system. Requires a function that can take a URL and load it
        into a file storage system, returning an ID from which the file can be
        retrieved. Generally you would get this from something like:
        
            with toil.common.Toil(options) as toil_instance:
                plan.bake(lambda url: toil_instance.importFile(url))
        
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
            imported_id = import_function(primary_url)
            Logger.info("Imported primary {} as {}".format(primary_url,
                imported_id))
            self.primary_ids.append(imported_id)
            
        for alt_url in self.alt_fasta_urls:
            # Import all the alt FASTAs
            imported_id = import_function(alt_url)
            Logger.info("Imported alt {} as {}".format(alt_url,
                imported_id))
            self.alt_ids.append(imported_id)
            
        for chrom_name, urls in self.vcf_urls.iteritems():
            # Import all the VCFs
            
            for vcf_url in urls:
                imported_id = import_function(vcf_url)
                Logger.info("Imported VCF {} as {}".format(vcf_url, imported_id))
                self.vcf_ids[chrom_name].append(imported_id)
                
                if self.vcf_index_urls.has_key(vcf_url):
                    # There's also an index
                    vcf_index_url = self.vcf_index_urls[vcf_url]
                    index_id = import_function(vcf_index_url)
                    Logger.info("Imported index {} as {}".format(vcf_index_url,
                        index_id))
                    # Remember that this is the index for that VCF
                    self.index_ids[imported_id] = index_id
            
        for hal_url in self.hal_urls:
            # Import all the HALs
            imported_id = import_function(hal_url)
            Logger.info("Imported HAL {} as {}".format(hal_url,
                imported_id))
            self.hal_ids.append(imported_id)
            
        for base_vg_url in self.base_vg_urls:
            # Import all the base VGs
            imported_id = import_function(base_vg_url)
            Logger.info("Imported base VG {} as {}".format(base_vg_url,
                imported_id))
            self.base_vg_ids.append(imported_id)
    
    def set_chromosome_name(self, accession, name):
        """
        Set the chromosome name for a primary scaffold. Takes the
        accession.version string (like "LOL1234.1") and the chromosome name
        string (like "1" or "MT").
        """

        self.accession_to_name[accession] = name
        self.name_to_accession[name] = accession
        
    def set_alt_parent(self, alt_accession, parent_accession):
        """
        Takes the accession.version string (like "LOL1234.1") for an alt
        scaffold and a primary scaffold, and makes the primary scaffold the
        parent of the alt scaffold. This means that the alt scaffold will be
        aligned against and merged into the subgraph derived from the primary
        scaffold and ist other alt scaffolds.
        """

        self.alt_parents[alt_accession] = parent_accession
        
    def get_alt_parent(self, alt_accession):
        """
        Takes the accession.version string (like "LOL1234.1") for an alt
        scaffold, and returns the accession.version string of the primary
        scaffold it belongs to.
        """

        return self.alt_parents[alt_accession]
        
    def get_children(self, parent_accession):
        """
        Iterate through the list of accession.version strings of child alt
        contigs on the given parent, identifird by accession.version.
        """
        
        # We'll build a list of the children of this parent
        to_return = []
        
        # TODO: reverse this so we have efficient access in this direction.
        for child, parent in self.alt_parents.iteritems():
            if parent == parent_accession:
                # This is a child of that parent
                to_return.append(child)
        
        return to_return
        
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
            self.set_chromosome_name(accession, name)
            
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
            
                   
        
        
        














    
    
