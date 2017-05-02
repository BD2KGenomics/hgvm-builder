#hgvm-builder commandrunner.py: tool for running commands and pipelines

import logging
import subprocess
import os
import os.path

# Get a submodule-global logger
Logger = logging.getLogger("commandrunner")

class CommandRunner(object):
    """
    Class to run pipelines of commands. Will eventually support running through Docker.
    
    Based on DockerRunner from <https://github.com/BD2KGenomics/toil-vg/blob/master/src/toil_vg/vg_common.py>.
    
    TODO: Pull in toil-vg as a dependency
    """

    def call(self, args, work_dir = '.' , outfile = None, errfile = None,
             check_output = False, inputs=[]):
        """
        Call the given command or pipeline directly.
        """

        Logger.info("Run: {}".format(" | ".join(" ".join(x) for x in args)))

        # this is all that docker_call does with the inputs parameter:
        for filename in inputs:
            assert(os.path.isfile(os.path.join(work_dir, filename)))

        procs = []
        for i in range(len(args)):
            stdin = procs[i-1].stdout if i > 0 else None
            if i == len(args) - 1 and outfile is not None:
                stdout = outfile
            else:
                stdout = subprocess.PIPE
            procs.append(subprocess.Popen(args[i], stdout=stdout, stderr=errfile,
                                          stdin=stdin, cwd=work_dir))
            
        for p in procs[:-1]:
            p.stdout.close()

        output, errors = procs[-1].communicate()
        for i, proc in enumerate(procs):
            sts = proc.wait()
            if sts != 0:            
                raise Exception("Command {} returned with non-zero exit status {}".format(
                    " ".join(args[i]), sts))

        if check_output:
            return output


