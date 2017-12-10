import os
import sys
import json
import multiprocessing_on_dill as multiprocessing

import fieldpathogenomics
from fieldpathogenomics.pipelines.Callset import VCFtoHDF5RefSnps
import fieldpathogenomics.utils as utils

from bioluigi.slurm import SlurmExecutableTask, SlurmTask
from bioluigi.utils import CheckTargetNonEmpty
from bioluigi.decorators import requires, inherits


import luigi
from luigi import LocalTarget

FILE_HASH = utils.file_hash(__file__)
PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]


@requires(VCFtoHDF5RefSnps)
class GetMSA(SlurmTask, CheckTargetNonEmpty):
    gene_cov = luigi.FloatParameter()
    gff = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 20000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, str(self.gene_cov), self.output_prefix + ".phy")),
                LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, str(self.gene_cov), self.output_prefix + ".sites"))]

    def work(self):
        import allel
        import h5py
        import Bio, Bio.AlignIO, Bio.Align
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from collections import defaultdict
        from fieldpathogenomics.utils import GeneSet, index_variants

        refsnps = h5py.File(self.input().path, mode='r')

        genotypes = allel.GenotypeDaskArray(refsnps['calldata']['GT'])
        variants = allel.VariantChunkedTable(refsnps['variants'])
        index = index_variants(refsnps['variants'], 7)
        samples = list(refsnps['samples'])
        gene_table = allel.FeatureTable.from_gff3(self.gff, attributes=['ID, Parent']).query("type == 'gene'")

        gs = GeneSet(gene_table, genotypes, variants, samples, index)

        msa = Bio.Align.MultipleSeqAlignment([SeqRecord(Seq(''), id=lib) for lib in samples])
        invariants = defaultdict(int)

        for cov, snps, inv in zip(gs.sample_coverages(), gs.MSA(), gs.invariants()):
            if cov.min() < self.gene_cov:
                continue

            msa += snps
            for b in 'ACGT':
                invariants[b] += inv.get(b, 0)

        with self.output()[0].open('w') as fphy:
            Bio.AlignIO.write(Bio.Align.MultipleSeqAlignment(msa), fphy, 'phylip-relaxed')
        with self.output()[1].open('w') as finv:
            finv.write('/'.join([str(invariants[b]) for b in 'ACGT']))


@requires(GetMSA)
class RAxML_ng(SlurmExecutableTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 8
        self.partition = "nbi-long,RG-Diane-Saunders"
        self.sbatch_args = '--constraint=intel'

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, str(self.gene_cov), 'ng', self.output_prefix + ".raxml.support"))

    def work_script(self):
        with self.input()[1].open('r') as f:
            asc = f.read().strip()
        return '''#!/bin/bash
               source raxmlng-0.5;

               rm -r {output_dir}_temp
               mkdir {output_dir}_temp
               cd {output_dir}_temp

               set -euo pipefail

               raxml-ng-mpi --all \
                            --redo \
                            --msa {input} \
                            --model GTR+G+ASC_STAM{{{asc}}} \
                            --prefix {prefix} \
                            --threads {n_cpu}

               mv {output_dir}_temp/* {output_dir}
               rmdir {output_dir}_temp
               '''.format(output_dir=os.path.split(self.output().path)[0],
                          n_cpu=self.n_cpu,
                          input=self.input()[0].path,
                          prefix=self.output_prefix,
                          asc=asc)


@inherits(RAxML_ng)
class CovWrapper(luigi.WrapperTask):
    gene_cov = None

    def requires(self):
        for c in [0.5, 0.6, 0.7, 0.8, 0.9]:
            yield self.clone(RAxML_ng, gene_cov=c)

    def output(self):
        return list(self.input())

if __name__ == '__main__':
    multiprocessing.set_start_method('forkserver')
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=os.path.basename(__file__))

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]

    luigi.run(['CovWrapper', '--output-prefix', name,
                           '--lib-list', json.dumps(lib_list),
                           '--star-genome', os.path.join(utils.reference_dir, 'genome'),
                           '--gff', os.path.join(utils.reference_dir, 'PST_genes_final.gff3'),
                           '--reference', os.path.join(utils.reference_dir, 'PST130_contigs.fasta'),
                           '--mask', os.path.join(utils.reference_dir, 'PST130_RNASeq_collapsed_exons.bed')] + sys.argv[2:])
