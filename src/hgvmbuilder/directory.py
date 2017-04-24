#hgvm-builder directory.py: Represent a bunch of Toil files with filenames

import os
import os.path
import logging
import urlparse
import shutil

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
        
    def __repr__(self):
        """
        Report ourselves fully.
        """
        
        return "Directory(" + repr(self.ids) + ")"
        
    def add(self, file_name, file_id):
        """
        Save the file with the given ID in the directory under the given name.
        No leading slashes should be used. "." and ".." are not supported.
        
        Can be chained.
        """
        
        self.ids[file_name] = file_id
        
        Logger.info("Saved {} as {}".format(file_id, file_name))
        
        return self
        
    def mount(self, prefix, other):
        """
        Add all the files from the other directory into this one with the given
        name prefix. A trailing slash will automatically be added if not
        present.
        
        Can be chained.
        """
        
        if not prefix.endswith("/"):
            # Add the trailing slash if not specified already.
            prefix = prefix + "/"
        
        for file_name, file_id in other.for_each_file():
            # Stick each file in under the prefix.
            self.add(prefix + file_name, file_id)
        
        return self
        
    def mount_all(self, dirs_by_prefix):
        """
        Mount all the directories in the given dict under their keys. Keys
        without a trailing "/" will have it added.
        
        Can be chained.
        
        """
        
        for prefix, directory in dirs_by_prefix.iteritems():
            self.mount(prefix, directory)
        
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
            
    def has_file(self, file_name):
        """
        Returns true if we have a file with the given name, and false otherwise.
        
        """
        
        return self.ids.has_key(file_name)
        
            
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
                
    def merge(self, *others):
        """
        Merge all the files in the given other directories into this one.
        
        Can be chained.
        """
        
        for other in others:
            for file_name, file_id in other.for_each_file():
                # Just loop through the files and add them
                self.add(file_name, file_id)
            
        return self
                
    def export_to(self, export_function, base_url):
        """
        Export all the files in the Cirectory under the given base URL, using
        the given exporting function that takes an ID and a target URL. Must be
        run on the Toil master.
        
        Can be chained.
        """
        
        Logger.info("Export directory to {}".format(base_url))
        
        # We need to parse the URL so we can do special handling of file URLs
        parsing = urlparse.urlparse(base_url)
        
        for file_name, file_id in self.for_each_file():
            # Export every file under the given base URL.
            
            if parsing.scheme == "file":
                # Hack to fix Toil not creating directories needed to hold the
                # files it exports.
                
                # Work out what directory needs to hold this file
                parent_directory = os.path.dirname(os.path.join(parsing.path,
                    file_name))
                try:
                    # Then make sure it exists
                    os.makedirs(parent_directory)
                except OSError:
                    pass
                          
            # Now that we know the directory exists if it's needed, export the
            # file.
            export_function(file_id, "{}/{}".format(base_url,
                file_name))
                
        return self
        
    @staticmethod
    def import_from(import_function, base_url, file_names=None):
        """
        Import the files and directories under the given file URL into Toil
        using the given URL-to-file-ID uploading function, and return a
        Directory representing them.
        
        For non-file URLs, file_names must be specified. Every file name in the
        collection will be imported.
        
        For file URLs, the directory is always walked and imported in its
        entirety.
        
        TODO: accept manifests that list files
        
        """
        
        Logger.info("Import directory from {}".format(base_url))
        
        # Make an empty Directory
        directory = Directory()
        
        # We need to parse the URL so we can do special handling of file URLs
        parsing = urlparse.urlparse(base_url)
        
        if parsing.scheme == "file":
            # Work out what directory to look in
            base_path = os.path.abspath(parsing.path)
            
            for dir_path, dir_names, file_names in os.walk(base_path):
                # For every directory we find
                for file_name in file_names:
                    # For every file in that directory
                    
                    # Get its full absolute path
                    full_path = os.path.join(base_path, dir_path, file_name)
                    
                    # Import it
                    file_id = import_function("file:" + full_path)
                    
                    # Register it in the directory we're building, using the path
                    # relative from the path for the base directory.
                    directory.add(os.path.relpath(os.path.join(dir_path,
                        file_name), base_path), file_id)
            
        elif file_names is not None:
            # Just try all the file names
            
            for file_name in file_names:
                # Import every file that was named
                directory.add(file_name,
                    import_function(base_url + "/" + file_name))
        else:
            raise RuntimeError("Need file name hints for non-file URLs!")
        
        # Return the completed directory
        return directory
            
    
    def dump(self, file_store, base_path):
        """
        Dump the directory from the given file store to the given filesystem
        path, so that the file store can't clean it up and it actually sticks
        around on disk.
        
        The target path should be absolute, because the job's working directory
        is probably in some kind of sandbox.
        
        Can be chained.
        """
        
        Logger.info("Dumping directory to {}".format(base_path))
        
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
            
            # Actually download the file to a temp directory.
            temp_name = file_store.readGlobalFile(file_id)
            
            # Copy it to the right place, ignoring permissions.
            shutil.copyfile(temp_name, local_name)
            
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
            
            # Actually download the file to the given directory. TODO: the
            # file_store may delete it later!
            file_store.readGlobalFile(file_id, local_name)
            
        return self
        
            
    
