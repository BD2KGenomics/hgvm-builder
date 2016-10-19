#hgvm-builder vcfrewriter.py: Allow VCFs to be rewritten to new coordinates

from collections import OrderedDict
import logging
import re

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
        
    def fix_stream(self, input_stream):
        """
        Scan the given (possibly compressed) input stream for fixable VCF
        errors, and fix them. Yields all lines as output.
        
        Currently only fixes missing description quotes in INFO lines, since the
        1000 Genomes VCFs have instances of that error.
        """
        
        # We handle our decompression ourselves, transparently, based on stream
        # content
        decompressed_stream = TransparentUnzip(input_stream)
        
        # We're looking for broken INFO lines that don't have leading quotes on
        # descriptions.
        info_desc_missing_leading_quote = re.compile(
            r'''(##INFO=<.*Description=)([^"][^"]*".*>.*)''',
            re.DOTALL)
        
        for line in decompressed_stream:
            # Check each line for this error
            match = info_desc_missing_leading_quote.match(line)
            
            if match is not None:
                # Insert the missing quote
                line = match.group(1) + "\"" + match.group(2)
                
            # Give out all the fixed-up lines
            yield line
        
        
    def rewrite_stream(self, input_stream, output_stream):
        """
        Rewrite the VCF data from the given possibly-compressed VCF data stream
        to the given output stream, in uncompressed format.
        
        Fixes fixable errors and renames contigs.
        """
        
        # Fix (and maybe decompress) the input stream
        fixed_stream = self.fix_stream(input_stream)
        
        # Make a VCF reader
        reader = vcf.Reader(fsock=fixed_stream,
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
        
