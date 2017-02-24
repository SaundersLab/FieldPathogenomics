import unittest
import luigi
import os

from fieldpathogenomics.luigi.notebook import NotebookTask

test_dir = os.path.split(__file__)[0]

class TestNB(NotebookTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100
        self.vars_dict = {"test_var": 'Hello World!'}


    def output(self):
        return luigi.LocalTarget(os.path.join(test_dir, 'scratch', "SLURM_TestOk.txt"))

    def work(self):
        print("WOOOOOO HELLO")
        with self.output().open('w') as f:
            f.write("OK\n")



class TestNotebookTask(unittest.TestCase):
    def test_Ok(self):
        task = TestNB(notebook='notebooks/test.ipynb')
        luigi.build([task], local_scheduler=True)
