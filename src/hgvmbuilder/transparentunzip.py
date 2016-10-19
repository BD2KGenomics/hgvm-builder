#hgvm-builder transparentunzip.py: Ungzip streams transparently if needed

import zlib

class TransparentUnzip(object):
    """
    A class that represents a transparently-un-gzipping stream.
    
    Only decompresses if necessary.
    """
    
    def __init__(self, stream):
        """
        Encapsulate the given stream in an unzipper.
        """
        
        # Keep the stream
        self.stream = stream
        
        # Make a decompressor with the accept a gzip header flag.
        # See <http://stackoverflow.com/a/22311297/402891>.
        # Will get None'd out after we're done decompressing and it is flushed.
        self.decompressor = zlib.decompressobj(zlib.MAX_WBITS + 16)
        
        # Are we really compressed? This gets set to True if we successfully
        # read the first block, and False if we failed to read the first block.
        # If it's False, there's no need to try to decompress any more blocks.
        self.header_read = None
        
        # We need to do lines ourselves, so we need to keep a buffer
        self.line_buffer = ""
        
        # Track out throughput
        self.compressed_bytes = 0
        self.uncompressed_bytes = 0
        
    def __iter__(self):
        """
        Get an iterator over lines in the stream. We are our own iterator.
        """
        return self
        
    def next(self):
        """
        Function as an iterator over lines.
        """
        
        # Read a line
        line = self.readline()
        
        if line == "":
            # No trailing newline only happens at the end of the file
            raise StopIteration
            
        # If we have any data, return it
        return line
        
    def readline(self, max_bytes=None):
        """
        Return the next line with trailing "/n", or "" if there is no next line.
        
        Won't return more than max_bytes bytes.
        
        Doesn't support the size limiting argument that file objects generally
        support.
        """
        
        # See if we have a line to spit out
        newline_index = self.line_buffer.find("\n")
        
        if newline_index == -1:
            # No line is in the buffer
            # Go get more data from the stream (in 16 k blocks)
            compressed = self.stream.read(16 * 2 ** 10)
            
            # Track the amount read
            self.compressed_bytes += len(compressed)
            
            if compressed == "":
                if self.decompressor is not None and self.header_read:
                    # No more input data, but we did sucessfully start
                    # decompressing. Flush out the output data.
                    self.line_buffer += self.decompressor.flush()
                    self.decompressor = None
                    # Spit out a line from that
                    return self.readline()
                
                # We didn't find a newline, and there's no more data, and
                # nothing to flush.
                
                if len(self.line_buffer) < max_bytes:
                    # We can fill the max requested bytes and have some left
                    # over
                    
                    # Pull out the max bytes wanted as the "line"
                    line = self.line_buffer[0:max_bytes]
                    
                    # Make the rest be in the buffer still
                    self.line_buffer = self.line_buffer[max_bytes:]
                    
                else:
                    # They want our whole buffer
                    # Take what we have (with no trailing \n added)
                    line = self.line_buffer
                    # Clear the buffer so next time we return "" as desired.
                    self.line_buffer = ""
                
                # Track the uncompressed bytes
                self.uncompressed_bytes += len(line)
                
                return line
                
            # Otherwise we found more data
            
            if self.header_read is None:
                # Try decompressing the first block
                try:
                    decompressed = self.decompressor.decompress(compressed)
                    # We sucessfully found the headers we needed
                    self.header_read = True
                except zlib.error:
                    # Just skip decompressing; it's probably not actually
                    # compressed.
                    decompressed = compressed
                    # We looked and didn't find a valid compressed header
                    self.header_read = False
                    
            else:
                # We know if we should be compressed or not
                if self.header_read:
                    # We do need to decompress
                    decompressed = self.decompressor.decompress(compressed)
                else:
                    # We don't need to decompress at all
                    decompressed = compressed
            
            # Stick it in the buffer
            self.line_buffer += decompressed
            
            # Try again. TODO: may stack overflow if there aren't newlines ever
            # and nobody gave us a max line length.
            return self.readline(max_bytes)
        
        else:
            
            if max_bytes is not None:
                # Adjust the newline index to return no more than max_bytes
                # bytes. Since it's an included character in the line, we have
                # to take 1 off of max_bytes.
                newline_index = min(newline_index, max_bytes - 1)
        
            # We have a line. Grab it.
            line = self.line_buffer[0:newline_index + 1]
            
            # Pop it off
            self.line_buffer = self.line_buffer[newline_index + 1:]
            
            # Track the uncompressed bytes
            self.uncompressed_bytes += len(line)
            
            # Return it
            return line
            
    def start_stats(self):
        """
        Reset byte stat tracking
        """
        
        self.compressed_bytes = 0
        self.uncompressed_bytes = 0
        
    def get_stats(self):
        """
        Return the number of compressed, uncompressed bytes processed.
        """
        
        return self.compressed_bytes, self.uncompressed_bytes        
