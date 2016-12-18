import unittest
import luigi
import os
import subprocess

import fieldpathogenomics.pipelines.Library as Library

test_dir = os.path.split(__file__)[0]


class TestPerLibPipeline(unittest.TestCase):
    '''Does an end to end test of the pipeline using
       test files of ~100 reads mapping to a 2000bp gene region
      PST130_9996:2000-4000
    '''

    def setUp(self):
        self.params = {
            'star_genome': os.path.join(test_dir, 'data', 'test_genome'),
            'reference': os.path.join(test_dir, 'data', 'test_reference.fasta'),
            'base_dir': os.path.join(test_dir, 'output'),
            'scratch_dir': os.path.join(test_dir, 'scratch'),
            'read_dir': os.path.join(test_dir, 'data'),
            'library': 'test',

        }
        self.make_star_index()

    def make_star_index(self):
        # Is big so generate as needed
        os.makedirs(self.params['star_genome'], exist_ok=True)
        cmd = "source star-2.5.2a; STAR --runMode genomeGenerate  --genomeFastaFiles {ref} --genomeDir {genome} --genomeSAindexNbases 5".format(
            ref=self.params['reference'], genome=self.params['star_genome']
        )
        subprocess.call(cmd, shell=True)

    def test_locally(self):
        plp = Library.PerLibPipeline(run_locally=True, **self.params)
        luigi.build([plp], local_scheduler=True)
