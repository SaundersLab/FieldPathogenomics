import unittest
import subprocess

from fieldpathogenomics.luigi.uv import UVExecutableTask


class TestOk(UVExecutableTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100

    def work_script(self):
        return '''#!/bin/bash
                  set -euo pipefail
                  echo "OK"'''

    def on_success(self):
        # Hook callback to capture output
        self.caught_err = self._fetch_task_failures()
        super().on_success()


class TestFail(UVExecutableTask):
    def __init__(self):
        super().__init__()
        self.n_cpu = 1
        self.mem = 100

    def work_script(self):
        return '''#!/bin/bash
                  set -euo pipefail
                  echo FAIL
                  exit 1'''


class TestUV(unittest.TestCase):
    def test_Ok(self):
        task = TestOk()
        task.run()
        self.assertEqual(task.caught_err, "OK")

    def test_Fail(self):
        task = TestFail()
        with self.assertRaises(subprocess.CalledProcessError):
            task.run()


if __name__ == '__main__':
    unittest.main()
