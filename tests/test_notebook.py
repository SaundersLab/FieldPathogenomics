import unittest
import luigi
import os

from fieldpathogenomics.luigi.notebook import NotebookTask

test_dir = os.path.split(__file__)[0]


class TestNB(NotebookTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, vars_dict={"test": 'Hello World!'}, **kwargs)
        self.n_cpu = 1
        self.mem = 100

    def output(self):
        return luigi.LocalTarget(os.path.join(test_dir, 'scratch', "test.ipynb"))


class TestNotebookTask(unittest.TestCase):
    def test_Ok(self):
        task = TestNB(notebook='/usr/users/ga004/buntingd/FP_dev/dev/src/fieldpathogenomics/tests/notebooks/test.ipynb')
        luigi.build([task], local_scheduler=True)
