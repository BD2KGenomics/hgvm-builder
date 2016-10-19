#hgvm-builder vcfrewriter.py: Allow VCFs to be rewritten to new coordinates

from collections import OrderedDict
import logging

import vcf

from .transparentunzip import TransparentUnzip

Logger = logging.getLogger("plan")

class VcfRewriter(object):
    """
    Provides a class that rewrites VCFs according to a contig name translation
    table.
    """
    
    def __init__(self, translation):
        """
        Create a new VcfRewriter. Translation must be a dict from old contig
        name to new contig name; contigs that are not mentioned will produce an
        error if they are encountered in a VCF.
        """
        
        self.translation = translation
        
    def rewrite_stream(self, input_stream, output_stream):
        """
        Rewrite the VCF data from the given possibly-compressed VCF data stream
        to the given output stream, in uncompressed format.
        """
        
        # We handle our decompression ourselves, transparently, based on stream
        # content
        decompressed_stream = TransparentUnzip(input_stream)
        
        # Make a VCF reader
        reader = vcf.Reader(fsock=decompressed_stream,
            strict_whitespace=True)
            
        self.rewrite_reader(reader, output_stream)
        
    def rewrite_reader(self, reader, output_stream):
        """
        Rewrite the VCF data from the given pyvcf Reader, writing to the given
        output stream.
        """
        
        # Rewrite all the contigs
        new_contigs = OrderedDict()
        for name, metadata in reader.contigs.iteritems():
            new_contigs[self.translation[name]] = metadata
            
        # Patch up the reader so it can be a template
        reader.contigs = new_contigs
        
        # Make a writer to spit out the translated VCF
        writer = vcf.Writer(output_stream, reader)
        
        for record in reader:
            # For every record in the input, fix up its CHROM field and write it
            # to the output.
            record.CHROM = self.translation[record.CHROM]
            writer.write_record(record)
            
        # We're done writing to the output now, so close it up cleanly.
        writer.flush()
        writer.close()
        
