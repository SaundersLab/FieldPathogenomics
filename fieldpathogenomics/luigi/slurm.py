import os,sys, re
import subprocess
import time
import logging
import random
import tempfile

import luigi

alloc_log = logging.getLogger('alloc_log')
logger = logging.getLogger('luigi-interface')

class SlurmMixin(object):
    '''Mixin to support running Task on a SLURM cluster
   
     Parameters:

    - n_cpu: Number of CPUs to allocate for the Task. 
    - mem: Amount of memory to require MB
    - partition: slurm partition to submit to
    - run_locally: Run locally instead of on the cluster.
    
     '''
    
    n_cpu = luigi.IntParameter(default=1, significant=False)
    mem = luigi.IntParameter(default=1000, significant=False)
    partition = luigi.Parameter(default='tgac-medium', significant=False)
    run_locally = luigi.BoolParameter(significant=False, description="run locally instead of on the cluster")
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

    def _salloc(self):
        '''Request a job allocation from the scheduler, blocks until its ready then return the job id '''
        salloc = "salloc -N 1 -c {n_cpu} -n 1 --mem {total_mem} -p {partition} -J {job_name} --no-shell".format(
        n_cpu=self.n_cpu, partition=self.partition, total_mem=int(self.mem*self.n_cpu), job_name=self.job_name)
        
        comp = subprocess.run(salloc, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, check=True)
        grant_id = re.compile('salloc: Granted job allocation (\S+)')
        
        for line in comp.stderr.split('\n'):
            if grant_id.match(line) is not None:
                return (grant_id.match(line).groups()[0])
                
        raise Exception("Unable to create job allocation: " + comp.stderr)

    def _srun(self, launch, alloc):
        '''Run the task in launch in allocation alloc'''
        srun = "srun -n 1 --kill-on-bad-exit  --quit-on-interrupt --jobid {jobid} -c {n_cpu} --mem-per-cpu {mem}  -o {outfile} -e {errfile} {launch}".format(
        n_cpu=self.n_cpu, jobid=alloc, mem=self.mem, launch=launch, outfile=self.outfile, errfile=self.errfile )
        ret = subprocess.run(srun, shell=True, check=True)

    def _slaunch(self, launch):
        return "salloc --quiet -N 1 -c {n_cpu} -n 1 --mem {total_mem} -p {partition} -J {job_name}  srun  -n 1 -c {n_cpu} --mem-per-cpu {mem} {launch} > {outfile} 2> {errfile}".format(n_cpu=self.n_cpu,
         mem=self.mem, partition=self.partition, job_name=self.job_name, launch=launch, outfile=self.outfile, errfile=self.errfile )

    def _fetch_task_failures(self):
        '''Handles reading either the local CompletedProcess or the SLURM log files'''
        ret = ''
        if self.run_locally:
            if not self.completedprocess.stdout is None:
                ret+= self.completedprocess.stdout.replace("\n", "\nstdout " + self.task_id + ": ") 
            if not self.completedprocess.stderr is None:
                ret+= self.completedprocess.stderr.replace("\n", "\nstdout " + self.task_id + ": ")
        else:
            try:
                with open(self.errfile, 'r') as err:
                    ret +="\nSLURM err " + self.task_id + ": " + err.read().replace("\n", "\nSLURM err " + self.task_id + ": ") 
            except (FileNotFoundError,AttributeError):
                ret +="\nSLURM err " + self.task_id + ": " + "None"
            try: 
                with open(self.outfile, 'r') as out:
                    ret +="\nSLURM out " + self.task_id + ": " + out.read().replace("\n", "\nSLURM out " + self.task_id + ": ") 
            except (FileNotFoundError,AttributeError):
                ret +="\nSLURM out " + self.task_id + ": " + "None"
            
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
    
    def scancel(self):
        if self.alloc is not None:
            subprocess.run("scancel {0}".format(self.alloc), shell=True, check=False)
class SlurmExecutableTask(luigi.Task, SlurmMixin):

    """
        Run an external executable on SLURM
        Override ``work_script()`` to return a shell script file as a string to run

    """

    def __init__(self, *args, **kwargs):
        super(SlurmExecutableTask, self).__init__(*args, **kwargs)
        self.job_name = self.task_family

    def run(self):
        self._init_tmp()            
        self.launcher = os.path.join(self.tmp_dir, "launch.sh")
        
        with open(self.launcher, 'w') as l:
            l.write(self.work_script())
            
        if self.run_locally:
            self.completedprocess = subprocess.run(self.launcher , shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            self.completedprocess.check_returncode()
            
        else:
            self.outfile = os.path.join(self.tmp_dir, 'job.out')
            self.errfile = os.path.join(self.tmp_dir, 'job.err')
            
            self.alloc = None
            try:
                self.alloc = self._salloc()
                logger.info("SLURM: jobid={0}".format(self.alloc) )
                alloc_log.info(self.task_id + "\t" + str(self.alloc))
                
                self._srun(self.launcher, self.alloc)
                
            finally:
                # Always be sure to free the slurm allocation
                self.scancel()

    def work_script(self):
        """Override this an make it return the shell script to run"""
        pass
