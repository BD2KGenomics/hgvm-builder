#hgvm-builder directory.py: Represent a bunch of Toil files with filenames fro dumping

import os
import os.path
import logging

Logger = logging.getLogger("directory")

class Directory(object):
    """
    Represents a collection of Toil file IDs with assigned file names (which can
    have directories in them). Allows a bunch of file IDs to be passed around as
    a group, dumped to disk, exported from the Toil master, and generally
    treated as a virtual filesystem.
    
    Note that changes are only visible to the current job, or child jobs to whom
    the modified Directory object is passed. Note that Toil may not allow you to
    retrieve a file you just wrote.
    
    File IDs may be Toil promises, but exporting and downloading will not work
    until they are filled in.
    
    """
    
    def __init__(self, values={}):
        """
        Make a new empty Directory with no files (by default), or the files
        given in the specified dictionary.
        
        """
        
        # This holds Toil IDs by file name
        self.ids = values
        
    def add(self, file_name, file_id):
        """
        Save the file with the given ID in the directory under the given name.
        No leading slashes should be used. "." and ".." are not supported.
        
        Can be chained.
        """
        
        self.ids[file_name] = file_id
        
        Logger.info("Saved {} as {}".format(file_id, file_name))
        
        return self
        
    def get(self, file_name):
        """
        Return the Toil file ID associated with the given file name, or throw an
        error if no such file exists.
        """
        
        if self.ids.has_key(file_name):
            return self.ids[file_name]
        else:
            raise Exception("File {} not found.".format(file_name))
            
    def for_each_file(self, prefix=""):
        """
        Iterate over name, ID pairs in the Directory that start with the given
        prefix.
        """
        
        for file_name, file_id in self.ids.iteritems():
            # Just loop over all the files
            if file_name.startswith(prefix):
                # And if we find one that matches, yield it
                yield (file_name, file_id)
                
    def merge(self, other):
        """
        Merge all the files in the given other directory into this one.
        
        Can be chained.
        """
        
        for file_name, file_id in other.for_each_file():
            # Just loop through the files and add them
            self.add(file_name, file_id)
            
        return self
                
    def export(self, toil_instance, base_url):
        """
        Export all the files in the Cirectory under the given base URL, using
        the given Toil instance. Must bve run on the Toil master.
        
        Can be chained.
        """
        
        Logger.info("Export directory to {}".format(base_url))
        
        for file_name, file_id in self.for_each_file():
            # Export every file under the given base URL. TODO: catch when Toil
            # is using file URLs and make it make subdirectories instead of
            # crashing.
            toil_instance.exportFile(file_id, "{}/{}".format(base_url,
                file_name))
                
        return self
                
    def download(self, file_store, base_path):
        """
        Download all the contents of the Directory to the given local path,
        which may not yet exist, using the given Toil FileStore. Can be run on
        any node, but writes only to the node's local filesystem.
        
        Can be chained.
        """
        
        Logger.info("Download directory to {}".format(base_path))
        
        for file_name, file_id in self.for_each_file():
            # For every file
            
            # Compute its destination local filename
            local_name = os.path.join(base_path, file_name)
            
            # Work out the directory part
            local_dirname = os.path.dirname(local_name)
            
            if not os.path.exists(local_dirname):
                # If the file needs to go in a directory that doesn't exist
                try:
                    # Make it
                    os.makedirs(local_dirname)
                except:
                    # But don't worry if someone else made it first
                    pass
            
            # Actually download the file to the given directory. TODO: can this
            # result in a broken symlink if the joibstore goes away?
            file_store.readGlobalFile(file_id, local_name)
            
        return self
        
            
    
