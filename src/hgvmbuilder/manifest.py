#hgvm-builder manifest.py: Filestore-serializable metadata

"""
Contains a Manifest class that acts like a dict of metadata (holding anything
JSON can hold) and which knows how to load itself from/save itself to a Toil
filestore.
"""

import os
import os.path
import logging
import itertools
import sys
import json

Logger = logging.getLogger("toilpromise")

class Manifest(dict):
    """
    Holds key-value metadata and can read and write itself to/from a filestore
    as JSON.
    """
    
    def __init__(self, contents={}):
        """
        Make a new Manifest from the given dict.
        
        """
        
        self.update(contents)
            
    @staticmethod
    def load(filestore, file_id):
        """
        Load a manifest from the given file ID in the given file store.
        """
    
        with filestore.readGlobalFileStream(file_id) as stream:
            # Read the contained object and make a Manifest of it
            return Manifest(json.load(stream))
            
                
    def save(self, filestore):
        """
        Save the serialized Manifest to the given Toil filestore, and return the
        resulting file ID.
        """
        
        with filestore.writeGlobalFileStream() as (file_handle, file_id):
            # Pretty-print to the file
            json.dump(self, file_handle, indent=2)
            # Ensure a trailing newline
            file_handle.write("\n")
        
        return file_id
        
