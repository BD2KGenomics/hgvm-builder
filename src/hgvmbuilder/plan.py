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
    
    def add_primary_scaffold_fasta(self, url):
        """
        Adds the (possibly gzipped) FASTA at the given URL as a FASTA known to
        contain only primary scaffolds (the assembled chromosomes, things like
        chrM, and other top-level unplaced scaffolds).
        
        FASTA records should be named as accession.version.
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
        
    def add_primary_scaffold_names(self, url):
        """
        Adds the TSV at the given URL, in GRC chr2acc format (chromosome
        number/X/Y/MT and accession.version, with #-comments) to the database of
        chromosome to accession.version mappings. This database is then used to
        interpret VCFs.
        """
        
        raise NotImplementedError
        
    def add_variants(self, chromosome_name, url): 
        """
        Adds the gzipped VCF at the given URL (with associated index at
        <url>.tbi) as variants to be applied to the given chromosome number or
        X/Y/MT name. Only one VCF is allowed per chromosome.
        """
        
        raise NotImplementedError
        
    def add_sample(self, fastq_url):
        """
        Adds the given high coverage sample FASTQ as a sample to be aligned to
        and used to augment the graph.
        
        TODO: Will be implemented after VCFs and assembled contigs work
        """
        
        raise NotImplementedError














    
    
