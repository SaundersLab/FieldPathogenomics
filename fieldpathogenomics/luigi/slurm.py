import os
import re
import subprocess
import luigi
import dill
import sys

from fieldpathogenomics.luigi.cluster import ClusterBase

import logging
alloc_log = logging.getLogger('alloc_log')
logger = logging.getLogger('luigi-interface')


class SlurmMixin(ClusterBase):
    '''Mixin to support running Task on a SLURM cluster

     Parameters:

    - n_cpu: Number of CPUs to allocate for the Task.
    - mem: Amount of memory to require MB
    - partition: slurm partition to submit to
    - run_locally: Run locally instead of on the cluster.

     '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _salloc(self):
        '''Request a job allocation from the scheduler, blocks until its ready then return the job id '''
        salloc = "salloc -N 1 -c {n_cpu} -n 1 --mem {total_mem} -p {partition} -J {job_name} --no-shell".format(
            n_cpu=self.n_cpu, partition=self.partition, total_mem=int(self.mem * self.n_cpu), job_name=self.job_name)

        comp = subprocess.run(salloc, shell=True, stderr=subprocess.PIPE,
                              stdout=subprocess.PIPE, universal_newlines=True, check=True)
        grant_id = re.compile('salloc: Granted job allocation (\S+)')

        for line in comp.stderr.split('\n'):
            if grant_id.match(line) is not None:
                return (grant_id.match(line).groups()[0])

        raise Exception("Unable to create job allocation: " + comp.stderr)

    def _srun(self, launch, alloc):
        '''Run the task in launch in allocation alloc'''
        srun = "srun -n 1 --kill-on-bad-exit  --quit-on-interrupt --jobid {jobid} -c {n_cpu} --mem-per-cpu {mem}  -o {outfile} -e {errfile} {launch}".format(
            n_cpu=self.n_cpu, jobid=alloc, mem=self.mem, launch=launch, outfile=self.outfile, errfile=self.errfile)
        subprocess.run(srun, shell=True, check=True)

    def _slaunch(self, launch):
        return "salloc --quiet -N 1 -c {n_cpu} -n 1 --mem {total_mem} -p {partition} -J {job_name}  srun  -n 1 -c {n_cpu} --mem-per-cpu {mem} {launch} > {outfile} 2> {errfile}".format(
            n_cpu=self.n_cpu, mem=self.mem, partition=self.partition, job_name=self.job_name, launch=launch, outfile=self.outfile, errfile=self.errfile)

    def scancel(self):
        if self.alloc is not None:
            subprocess.run("scancel {0}".format(self.alloc), shell=True, check=False)


class SlurmExecutableTask(luigi.Task, SlurmMixin):

    """
        Run an external executable on SLURM
        Override ``work_script()`` to return a shell script file as a string to run

    """
    n_cpu = luigi.IntParameter(default=1, significant=False)
    mem = luigi.IntParameter(default=1000, significant=False)
    partition = luigi.Parameter(default='tgac-medium', significant=False)
    run_locally = luigi.BoolParameter(
        significant=False, description="run locally instead of on the cluster")
    rm_tmp = luigi.BoolParameter(default=True, significant=False)

    def __init__(self, *args, **kwargs):
        super(SlurmExecutableTask, self).__init__(*args, **kwargs)
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
        pass


class SlurmTask(SlurmExecutableTask):

    def work_script(self):
        python = os.path.join((os.environ['VIRTUAL_ENV']), 'bin', 'activate')
        return '''#!/bin/bash
                  source {python}
                  set -euo pipefail
                  python -m fieldpathogenomics.luigi.task_runner {task}
                  '''.format(python=python,
                             task=self.job_file)

    def _dump(self):
        """Dump instance to file."""
        with self.no_unpicklable_properties():
            self.job_file = os.path.join(self.tmp_dir, 'job-instance.pickle')
            if self.__module__ == '__main__':
                d = dill.dumps(self)
                module_name = os.path.basename(sys.argv[0]).rsplit('.', 1)[0]
                d = d.replace(b'(c__main__', b"(c" + module_name.encode())
                open(self.job_file, "wb").write(d)
            else:
                dill.dump(self, open(self.job_file, "wb"))

    def run(self):
        # Bit of a hack, _init_tmp() also gets called again inside super().run()
        self._init_tmp()
        self._dump()
        super().run()
