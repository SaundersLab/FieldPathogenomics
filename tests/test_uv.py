import unittest
import subprocess
import luigi
import os

from fieldpathogenomics.luigi.uv import UVExecutableTask

test_dir = os.path.split(__file__)[0]


class TestOk(UVExecutableTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100
        self.host = 'uv2k2'

    def output(self):
        return luigi.LocalTarget(os.path.join(test_dir, 'scratch', "UV_TestOk.txt"))

    def work_script(self):
        return '''#!/bin/bash
                  set -euo pipefail
                  echo OK > {output}
                  '''.format(output=self.output().path)

    def on_success(self):
        # Hook callback to capture output
        self.caught_err = self._fetch_task_failures()
        super().on_success()

    def on_failure(self, exception):
        # Hook callback to capture output
        self.caught_err = self._fetch_task_failures()
        super().on_failure(exception)


class TestFail(UVExecutableTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100
        self.host = 'uv2k2'

    def work_script(self):
        return '''#!/bin/bash
                  set -euo pipefail
                  echo FAIL
                  exit 1'''


class TestUV(unittest.TestCase):
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
