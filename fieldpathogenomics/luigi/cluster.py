import os
import subprocess
import tempfile

import logging
logger = logging.getLogger('luigi-interface')


class ClusterBase(object):
    '''
    '''

    def _init_tmp(self):
        # Set up temp folder in shared directory
        base_tmp_dir = tempfile.gettempdir()
        self.tmp_dir = os.path.join(base_tmp_dir, self.task_id)
        logger.info("Tmp dir: %s", self.tmp_dir)
        os.makedirs(self.tmp_dir, exist_ok=True)

    def clear_tmp(self):
        try:
            if (os.path.exists(self.tmp_dir) and self.rm_tmp):
                logger.debug('Removing temporary directory %s' % self.tmp_dir)
                subprocess.call(["rm", "-rf", self.tmp_dir])
        except:
            pass

    def _read_std(self):
        stdout, stderr = None, None
        if self.run_locally:
            try:
                if self.completedprocess.stdout is not None:
                    stdout = self.completedprocess.stdout
                if self.completedprocess.stderr is not None:
                    stderr = self.completedprocess.stderr
            except AttributeError:
                pass
        else:
            try:
                with open(self.errfile, 'r') as err:
                    stderr = err.read()
            except (FileNotFoundError, AttributeError):
                stderr = None

            try:
                with open(self.outfile, 'r') as out:
                    stdout = out.read()
            except (FileNotFoundError, AttributeError):
                stdout = None
        return stderr, stdout

    @property
    def stdout(self):
        _, stdout = self._read_std()
        if stdout is not None:
            self._stdout = stdout
        try:
            return self._stdout
        except AttributeError:
            self._stdout = None
            return self._stdout

    @property
    def stderr(self):
        stderr, _ = self._read_std()
        if stderr is not None:
            self._stderr = stderr
        try:
            return self._stderr
        except AttributeError:
            self._stderr = None
            return self._stderr

    def format_log(self):
        '''Handles reading either the local CompletedProcess or the SLURM log files'''
        ret = ''
        if self.run_locally:
            if self.stderr is not None:
                ret += self.stderr.replace("\n", "\nstderr " + self.task_id + ": ")

            if self.stdout is not None:
                ret += self.stdout.replace("\n", "\nstdout " + self.task_id + ": ")

        else:
            ret += "\nSLURM err " + self.task_id + ": "
            if self.stderr is not None:
                ret += self.stderr.replace("\n", "\nSLURM err " + self.task_id + ": ")

            ret += "\nSLURM out " + self.task_id + ": "
            if self.stdout is not None:
                ret += self.stdout.replace("\n", "\nSLURM out " + self.task_id + ": ")

        return ret
