import os
import subprocess
import luigi

from fieldpathogenomics.luigi.cluster import ClusterBase

import logging
logger = logging.getLogger('luigi-interface')
alloc_log = logging.getLogger('uv-alloc')


class PBSMixin(ClusterBase):
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

    def _qsub(self, launch):
        command = "qsub -l select=1:mem={mem}MB:ncpus={n_cpu} -q Test -W block=true -o {outfile} -e {errfile} -N {job_name} {launch} ".format(
            n_cpu=self.n_cpu, mem=self.mem, job_name=self.job_name, launch=launch, outfile=self.outfile, errfile=self.errfile)
        p = subprocess.run(self._ssh(command), shell=True, check=True, stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE, universal_newlines=True)
        # 208067.UV00000010-P002
        return p.stdout.strip()

    def _ssh(self, command):
        return "ssh {host} '{command}' ".format(host=self.host, command=command)

    def qdel(self):
        if self.alloc is not None:
            subprocess.run(self._ssh("qdel {alloc}".format(alloc=self.alloc)), shell=True, check=False)


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
                self.alloc = self._qsub(self.launcher)
                logger.info("UV: jobid={0}".format(self.alloc))
                alloc_log.info(self.task_id + "\t" + str(self.alloc))

            finally:
                # Always be sure to free the slurm allocation
                self.qdel()

    def on_failure(self, exception):
        err = self.format_log()
        self.clear_tmp()
        logger.info(err)
        super_retval = super().on_failure(exception)
        ret = err if super_retval is None else err + "\n" + super_retval
        return ret

    def on_success(self):
        err = self.format_log()
        self.clear_tmp()
        logger.info(err)
        super_retval = super().on_success()
        ret = err if super_retval is None else err + "\n" + super_retval
        return ret

    def work_script(self):
        """Override this an make it return the shell script to run"""
