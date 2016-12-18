
import luigi
from luigi import Target, LocalTarget
from luigi.util import task_wraps, inherits, requires
from luigi.task import getpaths
import six
import os
import math
import copy


def indextarget(struct, idx):
    """
    Maps all Targets in a structured output to an indexed temporary file
    """
    if isinstance(struct, Target):
        base, *ext = struct.path.split('.', maxsplit=1)
        if len(ext) > 0:
            return LocalTarget(base + "_" + str(idx) + "." + ext[0])
        else:
            return LocalTarget(base + "_" + str(idx))
#    elif isinstance(struct, list):
#        indextarget(struct[0], idx)
    else:
        raise NotImplemented()


class ScatterGather():
    '''Decorator to transparently add Scatter-Gather parallelism to a Luigi task
    :param scatterTask must inherit and implement a run() method which maps
           a single input() file to an array of output() files
    :param scatterTask must inherit and implement a run() method which maps
           an array of input() files to a single output() file
    :param N the number of parts to scatter into

    Example
    =======

    class scatter(luigi.Task):
        def run(self):
            with self.input().open() as fin:
                inp = fin.readlines()
            perfile = math.ceil(len(inp)/len(self.output()))        
            for i,out in enumerate(self.output()):
                with out.open('w') as fout:
                    fout.writelines(inp[i*perfile:(i+1)*perfile])


    class gather(luigi.Task):
        def run(self):
            with self.output().open('w') as fout:
                for i in self.input():
                    with i.open('r') as fin:
                        fout.write(fin.read())


    @ScatterGather(scatter, gather, 10) 
    class ToBeScattered(luigi.Task):
        def run(self):
            with self.input().open('r') as fin:
                with self.output().open('w') as fout:
                    for l in fin:
                        fout.write("Done! " + l)

    '''

    def __init__(self, scatterTask, gatherTask, N):
        self.scatterTask = scatterTask
        self.gatherTask = gatherTask
        self.N = N

    def metaProgScatter(self, scattertask):
        meta_self = self

        @inherits(self.workTask)
        class Scatter(scattertask):

            def requires(self):
                wt_req = meta_self.workTask.requires(self)
                return wt_req[0] if isinstance(wt_req, list) else wt_req

            def output(self):
                return [indextarget(meta_self.workTask.input(self), i) for i in range(meta_self.N)]

            def to_str_params(self, only_significant=False):
                sup = super().to_str_params(only_significant)
                extras = {'input': self.input().path, 'N': str(meta_self.N)}
                return dict(list(sup.items()) + list(extras.items()))

        return Scatter

    def metaProgWork(self, worktask):
        meta_self = self

        class Work(worktask):
            SG_index = luigi.IntParameter()

            @property
            def task_family(self):
                return worktask.task_family

            def requires(self):
                s = self.clone(meta_self.Scatter)
                return [s] + super().requires()[1:] if isinstance(super().requires(), list) else [s]

            def input(self):
                inp = super().input()
                scattered = inp[0][self.SG_index]
                return [inp[0][self.SG_index]] + inp[1:] if isinstance(super().requires(), list) else inp[0][self.SG_index]

            def output(self):
                return indextarget(super().output(), self.SG_index)
        return Work

    def metaProgGather(self, gathertask):
        meta_self = self

        @inherits(self.workTask)
        class Gather(gathertask):
            SG_index = None

            def requires(self):
                return [self.clone(meta_self.Work, SG_index=i) for i in range(meta_self.N)]

            def output(self):
                return meta_self.workTask.output(self)

            def to_str_params(self, only_significant=False):
                sup = super().to_str_params(only_significant)
                extras = {'N': str(meta_self.N), 'output': self.output().path}
                return dict(list(sup.items()) + list(extras.items()))

        return Gather

    def __call__(self, workTask):

        self.workTask = workTask

        self.Scatter = self.metaProgScatter(self.scatterTask)
        self.Work = self.metaProgWork(self.workTask)
        self.Gather = self.metaProgGather(self.gatherTask)

        return self.Gather
