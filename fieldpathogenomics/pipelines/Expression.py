import os
import sys
import json
import multiprocessing_on_dill as multiprocessing

import luigi
from luigi import LocalTarget

from bioluigi.slurm import SlurmExecutableTask, SlurmTask
from bioluigi.utils import CheckTargetNonEmpty
from bioluigi.decorators import requires, inherits

import fieldpathogenomics
import fieldpathogenomics.utils as utils
import fieldpathogenomics.pipelines.Library as Library


FILE_HASH = utils.file_hash(__file__)
VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]
PIPELINE = os.path.basename(__file__).split('.')[0]


class Transcriptome(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Pull out the spliced exons from the genomes'''
    gff = luigi.Parameter()
    reference = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(significant=False)


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'spliced_exons'))

    def work_script(self):
        return '''#!/bin/bash -e
                source gffread-0.9.8;
                source kallisto-0.43.0

                gffread {gff} -g {input} -w /dev/stdout | fold -w 60 > {output}.fasta
                kallisto index -i {output}.temp {output}.fasta

                mv {output}.temp {output}
                '''.format(gff=self.gff,
                           input=self.reference,
                           output=self.output().path)


@requires(index=Transcriptome, library=Library.Trimmomatic)
class KallistoQuant(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "nbi-short"

    def output(self):
        return {'tsv': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, 'abundance.tsv')),
                'h5': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, 'abundance.h5'))}

    def work_script(self):
        return '''#!/bin/bash -e
                source kallisto-0.43.0

                kallisto quant -t {n_cpu} -o {output_dir}_temp -b 100 -i {index} {R1} {R2}

                mv {output_dir}_temp/* {output_dir}
                rmdir {output_dir}_temp
          '''.format(output_dir=os.path.dirname(self.output()['tsv'].path),
                     index=self.input()['index'].path,
                     R1=self.input()['library'][0].path,
                     R2=self.input()['library'][1].path,
                     n_cpu=self.n_cpu)


@inherits(KallistoQuant)
class AggregateKallisto(CheckTargetNonEmpty, SlurmTask):
    library = None
    output_prefix = luigi.Parameter()
    lib_list = luigi.ListParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.partition = 'main'
        self.mem = 8000
        self.n_cpu = 1

    def requires(self):
        return {lib: self.clone(KallistoQuant, library=lib) for lib in self.lib_list}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, "kallisto_abundances.tsv"))

    def run(self):
        import pandas as pd
        import h5py
        import numpy as np
        from luigi.file import atomic_file

        d = {}
        for lib, p in self.input().items():
            with h5py.File(p['h5'].path, mode='r') as f:
                tpms = np.stack([f['bootstrap'][k][:] for k in f['bootstrap'].keys()])
            d[lib + '_tpm'] = tpms.mean(axis=0)
            d[lib + '_std'] = tpms.std(axis=0)

        with h5py.File(list(self.input().values())[0]['h5'].path, mode='r') as f:
            index = [x.decode() for x in f['aux']['ids']]

        df = pd.DataFrame(d, index = index)
        df.columns = pd.MultiIndex.from_tuples([tuple(c.split('_')) for c in df.columns])

        af = atomic_file(self.output().path)
        df.to_csv(af.tmp_path, sep='\t')
        af.move_to_final_destination()

# ----------------------------------------------------------------------- #


if __name__ == '__main__':
    multiprocessing.set_start_method('forkserver')
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=os.path.basename(__file__))

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]

    luigi.run(['AggregateKallisto', '--output-prefix', name,
                                    '--lib-list', json.dumps(lib_list),
                                    '--gff', os.path.join(utils.reference_dir, 'PST_genes_final.gff3'),
                                    '--reference', os.path.join(utils.reference_dir, 'PST130_contigs.fasta')] + sys.argv[2:])

