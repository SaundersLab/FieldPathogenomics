import unittest
import subprocess
import luigi
import os

from bioluigi.slurm import SlurmExecutableTask, SlurmTask

test_dir = os.path.split(__file__)[0]


class TestExecOk(SlurmExecutableTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100

    def output(self):
        return luigi.LocalTarget(os.path.join(test_dir, 'scratch', "SLURM_TestExecOk.txt"))

    def work_script(self):
        return '''#!/bin/bash
                  set -euo pipefail
                  echo OK > {output}
                  '''.format(output=self.output().path)


class TestExecFail(SlurmExecutableTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100

    def work_script(self):
        return '''#!/bin/bash
                  set -euo pipefail
                  echo FAIL
                  exit 1'''


class TestSlurmExecutableTask(unittest.TestCase):
    def test_Ok(self):
        task = TestExecOk()
        luigi.build([task], local_scheduler=True)
        with task.output().open() as f:
            self.assertEqual(f.readlines()[0], "OK\n")

    def test_Fail(self):
        task = TestExecFail()
        with self.assertRaises(subprocess.CalledProcessError):
            task.run()


class TestOk(SlurmTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100

    def output(self):
        return luigi.LocalTarget(os.path.join(test_dir, 'scratch', "SLURM_TestOk.txt"))

    def work(self):
        print("WOOOOOO HELLO")
        with self.output().open('w') as f:
            f.write("OK\n")


class TestFail(SlurmTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100

    def work(self):
        raise Exception()


class TestSlurmTask(unittest.TestCase):
    def test_Ok(self):
        task = TestOk()
        luigi.build([task], local_scheduler=True)
        with task.output().open() as f:
            self.assertEqual(f.readlines()[0], "OK\n")

    def test_Fail(self):
        task = TestFail()
        with self.assertRaises(subprocess.CalledProcessError):
            task.run()


if __name__ == '__main__':
    unittest.main()
