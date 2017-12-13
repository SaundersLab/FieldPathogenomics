import os
import sys
import json
import multiprocessing_on_dill as multiprocessing

import fieldpathogenomics
from fieldpathogenomics.pipelines.Tree import CovWrapper
from fieldpathogenomics.pipelines.Callset import HD5s
import fieldpathogenomics.utils as utils

from bioluigi.notebook import NotebookTask
from bioluigi.utils import CheckTargetNonEmpty
from bioluigi.decorators import requires, inherits

import luigi
from luigi import LocalTarget

FILE_HASH = utils.file_hash(__file__)
PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]


@requires(hd5s=HD5s, trees=CovWrapper)
class DAPC(NotebookTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mem = 16000
        self.n_cpu = 1
        self.partition = "nbi-medium"
        self.notebook = os.path.join(utils.notebooks, 'DAPC', 'DAPC.ipynb')
        self.vars_dict = {"SNP_HD5": self.input()['hd5s']['snps'].path,
                          "TREE_NWK": self.input()['trees'][0].path,
                          "MIN_COV": 0.5,
                          "CLUSTER_MIN": 2,
                          "N_CLUST": 4,
                          "CLUSTER_MAX": 10,
                          "MAX_LINKAGE": 0.95}
        logger.info(str(self.vars_dict))

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'DAPC.ipynb'))

# -----------------------------------------------------------------------------------

if __name__ == '__main__':
    multiprocessing.set_start_method('forkserver')
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=os.path.basename(__file__))

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]

    luigi.run(['DAPC', '--output-prefix', name,
                       '--lib-list', json.dumps(lib_list),
                       '--gff', os.path.join(utils.reference_dir, 'PST_genes_final.gff3'),
                       '--star-genome', os.path.join(utils.reference_dir, 'genome'),
                       '--reference', os.path.join(utils.reference_dir, 'PST130_contigs.fasta'),
                       '--mask', os.path.join(utils.reference_dir, 'PST130_RNASeq_collapsed_exons.bed')] + sys.argv[2:])
