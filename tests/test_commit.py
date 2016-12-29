import unittest
import luigi
import os

from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from fieldpathogenomics.luigi.commit import CommittedTarget, CommittedTask

test_dir = os.path.split(__file__)[0]


class TestTask(CommittedTask, SlurmExecutableTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100

    def output(self):
        return CommittedTarget(os.path.join(test_dir, 'scratch', "TestCommit.txt"))

    def run(self):
        with self.output().open('w') as f:
            f.write("Testing")


class TestCommit(unittest.TestCase):
    def test_Ok(self):
        task = TestTask()
        luigi.build([task], local_scheduler=True)


if __name__ == '__main__':
    unittest.main()
