import os
import subprocess
import tempfile
import luigi

import logging
logger = logging.getLogger('luigi-interface')
alloc_log = logging.getLogger('uv-alloc')


class PBSMixin(object):
    '''Mixin to support running Task on a SLURM cluster

     Parameters:

    - n_cpu: Number of CPUs to allocate for the Task.
    - mem: Amount of memory to require MB
    - run_locally: Run locally instead of on the cluster.

     '''

    n_cpu = luigi.IntParameter(default=1, significant=False)
    mem = luigi.IntParameter(default=1000, significant=False)
    host = luigi.Parameter(default='uv2k1', significant=False)
    run_locally = luigi.BoolParameter(
        significant=False, description="run locally instead of on the cluster")
    rm_tmp = luigi.BoolParameter(default=True, significant=False)

    def _init_tmp(self):
        # Set up temp folder in shared directory
        base_tmp_dir = tempfile.gettempdir()
        self.tmp_dir = os.path.join(base_tmp_dir, self.task_id)
        logger.info("Tmp dir: %s", self.tmp_dir)
        os.makedirs(self.tmp_dir)

    def clear_tmp(self):
        try:
            if (self.tmp_dir and os.path.exists(self.tmp_dir) and self.rm_tmp):
                logger.debug('Removing temporary directory %s' % self.tmp_dir)
                subprocess.call(["rm", "-rf", self.tmp_dir])
        except:
            pass

    def _qsub(self, launch):
        return "ssh {host} 'qsub -l select=1:mem={mem}MB:ncpus={n_cpu} -q Test -W block=true -keo -o {outfile} -e {errfile} -N {job_name} {launch} '".format(
            host=self.host, n_cpu=self.n_cpu, mem=self.mem, job_name=self.job_name, launch=launch, outfile=self.outfile, errfile=self.errfile)

    def _fetch_task_failures(self):
        '''Handles reading either the local CompletedProcess or the SLURM log files'''
        ret = ''
        if self.run_locally:
            if self.completedprocess.stdout is not None:
                ret += self.completedprocess.stdout.replace("\n", "\nstdout " + self.task_id + ": ")
            if self.completedprocess.stderr is not None:
                ret += self.completedprocess.stderr.replace("\n", "\nstdout " + self.task_id + ": ")
        else:
            try:
                with open(self.errfile, 'r') as err:
                    ret += "\nSLURM err " + self.task_id + ": " + \
                        err.read().replace("\n", "\nPBS err " + self.task_id + ": ")
            except (FileNotFoundError, AttributeError):
                ret += "\nSLURM err " + self.task_id + ": " + "None"
            try:
                with open(self.outfile, 'r') as out:
                    ret += "\nSLURM out " + self.task_id + ": " + \
                        out.read().replace("\n", "\nPBS out " + self.task_id + ": ")
            except (FileNotFoundError, AttributeError):
                ret += "\nSLURM out " + self.task_id + ": " + "None"

        return ret

    def on_failure(self, exception):
        err = self._fetch_task_failures()
        self.clear_tmp()
        logger.info(err)
        super_retval = super().on_failure(exception)
        if super_retval is not None:
            return err + "\n" + super_retval
        else:
            return err

    def on_success(self):
        err = self._fetch_task_failures()
        self.clear_tmp()
        logger.info(err)
        super_retval = super().on_success()
        if super_retval is not None:
            return err + "\n" + super_retval
        else:
            return err


class UVExecutableTask(luigi.Task, PBSMixin):

    """
        Run an external executable on SLURM
        Override ``work_script()`` to return a shell script file as a string to run

    """

    def __init__(self, *args, **kwargs):
        super(UVExecutableTask, self).__init__(*args, **kwargs)
        self.job_name = self.task_family

    def run(self):
        self._init_tmp()
        self.launcher = os.path.join(self.tmp_dir, "launch.sh")

        with open(self.launcher, 'w') as l:
            l.write(self.work_script())

        if self.run_locally:
            self.completedprocess = subprocess.run(
                self.launcher, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            self.completedprocess.check_returncode()

        else:
            self.outfile = os.path.join(self.tmp_dir, 'job.out')
            self.errfile = os.path.join(self.tmp_dir, 'job.err')

            self.alloc = None
            try:
                self.alloc = self._salloc()
                logger.info("SLURM: jobid={0}".format(self.alloc))
                alloc_log.info(self.task_id + "\t" + str(self.alloc))

                self._srun(self.launcher, self.alloc)

            finally:
                # Always be sure to free the slurm allocation
                self.scancel()

    def work_script(self):
        """Override this an make it return the shell script to run"""
