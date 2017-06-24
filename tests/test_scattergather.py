import unittest
import math
import os
import luigi
import luigi.mock
import glob

from bioluigi.scattergather import ScatterGather

test_dir = os.path.split(__file__)[0]


class scatter(luigi.Task):

    def run(self):
        with self.input().open('r') as fin:
            inp = fin.readlines()
        perfile = math.ceil(len(inp) / len(self.output()))
        for i, out in enumerate(self.output()):
            with out.open('w') as fout:
                fout.writelines(inp[i * perfile:(i + 1) * perfile])


class gather(luigi.Task):

    def run(self):
        with self.output().open('w') as fout:
            for i in self.input():
                with i.open('r') as fin:
                    fout.write(fin.read())


class testdata(luigi.ExternalTask):

    def output(self):
        t = luigi.LocalTarget(os.path.join(test_dir, "sgtest.txt"))
        with t.open('w') as f:
            f.write('\n'.join([str(x) for x in range(100)]))
        return t


class TestSimple(unittest.TestCase):

    def setUp(self):
        @ScatterGather(scatter, gather, 10)
        class Simple(luigi.Task):

            def requires(self):
                return testdata()

            def run(self):
                with self.input().open('r') as fin:
                    with self.output().open('w') as fout:
                        for l in fin:
                            fout.write("Done! " + l)

            def output(self):
                return luigi.LocalTarget(os.path.join(test_dir, "sgtest_out.txt"))

        self.task = Simple()

    def test_simple(self):

        self.assertTrue(luigi.build([self.task], local_scheduler=True))
        with self.task.output().open('r') as f:
            out = f.read()
        self.assertEqual(out, '\n'.join(["Done! " + str(x) for x in range(100)]))

    def tearDown(self):
        for f in glob.glob(os.path.join(test_dir, "sgtest*")):
            os.remove(f)


class TestMultiReqs(unittest.TestCase):

    def setUp(self):
        class OtherReq(luigi.ExternalTask):

            def output(self):
                t = luigi.LocalTarget(os.path.join(test_dir, "sgtest_other.txt"))
                with t.open('w') as f:
                    f.write('Done! ')
                return t

        @ScatterGather(scatter, gather, 10)
        class MultiReqs(luigi.Task):

            def requires(self):
                return [testdata(), OtherReq()]

            def run(self):

                with self.input()[0].open('r') as fin, self.input()[1].open('r') as f2in:
                    tag = f2in.read().strip()
                    with self.output().open('w') as fout:
                        for l in fin:
                            fout.write(tag + l)

            def output(self):
                return luigi.LocalTarget(os.path.join(test_dir, "sgtest_out.txt"))

        self.task = MultiReqs()

    def test_simple(self):

        self.assertTrue(luigi.build([self.task], local_scheduler=True))
        with self.task.output().open('r') as f:
            out = f.read()
        self.assertEqual(out, '\n'.join(["Done!" + str(x) for x in range(100)]))

    def tearDown(self):
        for f in glob.glob(os.path.join(test_dir, "sgtest*")):
            os.remove(f)

if __name__ == '__main__':
    unittest.main()
