#hgvm-builder ftputil.py: FTP client utility functions

import ftplib
import logging
import os
import os.path
import random
import select
import time
import urlparse

# Get a submodule-global logger
Logger = logging.getLogger("ftputil")

class FTPOrFilesystemConnection:
    """
    Represents a connection to an FTP server or to a directory on the local
    filesystem. When constructed from a root URL, allows listing of
    files/directories under relative paths from that URL, as well as getting a
    URL to a relative path.
    """
    
    def __init__(self, base_url):
        """
        Create a new FTPOrFilesystemConnection to an ftp:// or file:// URL.
        """
        
        # Save the base URL
        self.base_url = base_url
        
        # Make sure it has a trailing / for urljoin to work and put things
        # inside of here.
        if not self.base_url.endswith("/"):
            self.base_url = self.base_url + "/"
        

        # Split up the URL
        parts = urlparse.urlparse(base_url)
        
        if parts.scheme == "ftp":
            # If it's an FTP URL open an FTP connection
            self.connection, self.base_path = ftp_connect(base_url)
            
            if not self.base_path.endswith("/"):
                # We need the trailing slash for urljoin to work
                self.base_path = self.base_path + "/"
            
            Logger.debug("Connected to {} with base path {}".format(
                base_url, self.base_path))
        elif parts.scheme == "file":
            # If it's a File URL make a FakeFTP
            self.connection = FakeFTP(parts.path)
            # We don't have any base path within that; we pretend this is the
            # FTP root.
            self.base_path = ""
        else:
            # Otherwise explode. TODO: can we unify this with IOstores from my
            # old toillib? Or would that pull in too much Azure/AWS code...
            raise RuntimeError("Unsupported protocol {}".format(parts.scheme))        
        
        
    def list_children(self, path):
        """
        Yield all direct children of the given root-relative path, as root-
        relative paths. If the path refers to a file, there will be no children.
        """
        
        if len(path) > 0 and not path.endswith("/"):
            # We need a trailing slash after a real directory name for urljoin
            # to work later
            path = path + "/"
        
        # Strip leading slashes from the input path, so we always look inside
        # our base path.
        path = path.lstrip("/")
        
        # Construct the path to actually go to on the FTP server
        ftp_path = urlparse.urljoin(self.base_path, path)
        
        Logger.debug("Listing {}".format(ftp_path))
        
        for child in robust_nlst(self.connection, ftp_path):
            # For every child, build a root-relative URL
            yield urlparse.urljoin(path, child)
        
    def get_url(self, path):
        """
        Returns the full URL to the given root-relative path.
        """
        
        # Strip leading slashes from the input path, so we always look inside
        # our base path.
        path = path.lstrip("/")
        
        # We know our base URL has a trailing slash, so we can just urljoin onto
        # it.
        return urlparse.urljoin(self.base_url, path)

def backoff_times(retries=float("inf"), base_delay=300):
    """
    A generator that yields times for exponential back-off. Always yields 0
    first, and then if you have nonzero retries it yields exponentially but
    randomly increasing times in seconds to wait before trying again, stopping
    at the specified number of retries (and continuing forever by default).
    
    You have to do the error catching and sleeping yourself.
    
    Raises an exception if you use up all your backoff times.
    """
    
    # Don't wait at all before the first try
    yield 0
    
    # What retry are we on?
    try_number = 1
    
    # Make a delay that increases
    delay = float(base_delay) * 2
    
    while try_number <= retries:
        # Wait a random amount between 0 and 2^try_number * base_delay
        yield random.uniform(base_delay, delay)
        delay *= 2
        try_number += 1
        
    raise RuntimeError("Ran out of retries")

def ftp_connect(url, retries=float("inf")):
    """
    Connect to an FTP server and go to the specified directory with FTPlib.
    
    Return the ftplib connection and the path.
    
    On errors like timeouts, retry up to the given number of retries (infinite
    by default) with exponential back-off.
    
    """
    
    for delay in backoff_times(retries=retries):
        if delay > 0:
            # We have to wait before trying again
            Logger.info("Retry after {} seconds".format(delay))
            time.sleep(delay)
        try:
        
            ftp_info = urlparse.urlparse(url)
            assert(ftp_info.scheme == "ftp")
            
            # Connect to the server
            ftp = ftplib.FTP(ftp_info.netloc, ftp_info.username,
                ftp_info.password)
            
            # Log in
            ftp.login()
                
            # Go to the right directory
            ftp.cwd(ftp_info.path)
            
            return ftp, ftp_info.path
            
        except IOError as e:
            # Something went wrong doing the IO
            Logger.warning("Retry after FTP setup IO error: {}".format(e))
            
class FakeFTP:
    """
    This fakes an FTP connection on a normal directory. You can connect to it
    and use robust_nlst on it and iterate over the filesystem and FTP with the
    same API.
    """
    
    def __init__(self, root):
        """
        Make a new FakeFTP on the given root.
        """
        
        # Where is the root of our fake FTP server
        self.root = root
        
        # We need to be able to change directories
        self.relative_path = ""
        
        Logger.debug("Fake FTP root: {}".format(self.root))
        
    def nlst(self):
        """
        Return a list of all the files in the current directory.
        """
        Logger.debug("Fake FTP root: {} cwd: {}".format(self.root,
            self.relative_path))
        
        dir_path = os.path.join(self.root, self.relative_path)
        
        if os.path.isdir(dir_path):
            # Only directories have anything in them
            return os.listdir(dir_path)
        else:
            return []
        
    def cwd(self, path):
        """
        Change to the given directory relative to the root.
        """
        self.relative_path = path.lstrip("/")

def robust_nlst(ftp, path):
    """
    Occasionally ftplib (or some servers?) can get confused when the multiple
    sockets in play get out of sync. This wrapper attempts to fix up those
    situations.
    
    Given an absolute path relative to FTP root, return the basenames of all
    files and directories inside that path. If the path is actually a file,
    returns an empty list.
    
    """

    # We'll populate this with file/directory names. Start it balnk in case we
    # can't cwd because we're trying to cwd to a file.
    listing = []

    try:
        ftp.cwd(path)
    
        retries_remaining = 3
        
        while retries_remaining > 0:
            retries_remaining -= 1
            try:
                listing = ftp.nlst()
                break
            except ftplib.error_proto:
                # For some reason this built-in function sometimes doesn't work.
                # Try doing it ourselves.
                
                Logger.warning("FTPlib nlst on {} failed. "
                    "Clearing and retrying.".format(path))
                
                while True:
                    # Clear out any crap that may be coming over the ftp
                    # connection. Use a 2 second timeout
                    flags = select.select([ftp.sock], [], [], 2)
                    
                    if not flags[0]:
                        # Nothing came
                        break
                        
                    # Otherwise read some data and select again
                    got = ftp.sock.recv(1024)
                    Logger.warning("Extra data: {}".format(got))
                
                # Now manually do the listing
                listing = []
                try:
                    ftp.retrlines("NLST", lambda line: listing.append(line))
                    break
                except ftplib.error_reply:
                    Logger.warning("Manual NLST on {} failed. "
                        "Clearing and retrying.".format(path))
                
                    # The server... replied to a PASV wrong or something? Maybe
                    # we're out of sync. Try syncing up again.
                    while True:
                        # Clear out any crap that may be coming over the ftp
                        # connection. Use a 2 second timeout
                        flags = select.select([ftp.sock], [], [], 2)
                        
                        if not flags[0]:
                            # Nothing came
                            break
                            
                        # Otherwise read some data and select again
                        got = ftp.sock.recv(1024)
                        Logger.warning("Extra data: {}".format(got))
                            
    except ftplib.error_perm as e:
        error_code = int(e.args[0][:3])
        if error_code == 550:
            # We expect to do a lot of CWD-ing to files, raising this
            pass
        else:
            raise e                    
    
    # Return the list of files we managed to get.
    return listing
