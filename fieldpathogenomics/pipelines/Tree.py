import os
import sys
import json
import multiprocessing

import fieldpathogenomics
from fieldpathogenomics.pipelines.Callset import GetRefSNPs
import fieldpathogenomics.utils as utils

from bioluigi.slurm import SlurmExecutableTask, SlurmTask
from bioluigi.utils import CheckTargetNonEmpty
from bioluigi.decorators import requires, inherits


import luigi
from luigi import LocalTarget

FILE_HASH = utils.file_hash(__file__)
PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]


@requires(GetRefSNPs)
class ConvertToBCF(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Use bcftools view to convert the vcf to bcf, its worth doing this conversion
     as the bcf formatted file is much faster for separating the samples than the vcf'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 2
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.output_prefix + ".bcf.gz"))

    def work_script(self):
        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                bcftools view {input} -o {output}.temp -O b --threads 1

                mv {output}.temp {output}
                '''.format(input=self.input().path,
                           output=self.output().path)


@requires(ConvertToBCF)
class GetSingleSample(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Pull a single sample out of the joint BCF file and use picard to compute a BED file where there is missing data'''
    library = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 2
        self.partition = "nbi-short"

    def output(self):
        return[LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, 'single_sample', self.library + ".vcf.gz")),
               LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, 'single_sample', self.library + ".bed"))]

    def work_script(self):

        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                source vcftools-0.1.13;
                source jre-8u92
                source picardtools-2.1.1
                picard='{picard}'

                bcftools view {input} -o {vcf}.temp -O z -s {library} --exclude-uncalled --no-update
                zcat {vcf}.temp | sed 's/<NON_REF>/N/' | bgzip -c > {vcf}.temp1
                rm {vcf}.temp
                mv {vcf}.temp1 {vcf}
                tabix -p vcf {vcf}

                $picard IntervalListTools I={vcf} INVERT=true O=/dev/stdout | $picard IntervalListToBed I=/dev/stdin O={mask}.temp
                mv {mask}.temp {mask}
                '''.format(picard=utils.picard.format(mem=self.mem),
                           input=self.input().path,
                           vcf=self.output()[0].path,
                           mask=self.output()[1].path,
                           library=self.library)


@requires(GetSingleSample)
class BCFtoolsConsensus(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Apply the variants in the vcf file and mask with the missing data BED file.
       :param str consensus_type: can be H1 or H2 for applying one of the (pseudo) haplotypes or iupac-codes for ambiguous coding'''

    consensus_type = luigi.Parameter()  # {'H1' , 'H2', 'iupac-codes'}
    reference = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, 'single_sample', self.library + "_" + self.consensus_type + '.fasta'))

    def work_script(self):
        ct_flag = "--" + self.consensus_type if self.consensus_type == 'iupac-codes' else "-" + \
            self.consensus_type
        return '''#!/bin/bash -e

                source bcftools-1.3.1;
                bcftools consensus {vcf} -f {reference} -s {library} -m {mask} {consensus_type} > {output}.temp

                mv {output}.temp {output}

                source samtools-1.3;
                samtools faidx {output}

                '''.format(vcf=self.input()[0].path,
                           mask=self.input()[1].path,
                           reference=self.reference,
                           library=self.library,
                           output=self.output().path,
                           consensus_type=ct_flag)


@requires(BCFtoolsConsensus)
class GFFread(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Pull out the spliced exons from the genomes'''
    gff = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, 'single_sample', self.library + '_genes_' + self.consensus_type + '.fasta'))

    def work_script(self):
        return '''#!/bin/bash -e
                source gffread-0.9.8;

                gffread {gff} -g {input} -x /dev/stdout | fold -w 60 > {output}.temp

                mv {output}.temp {output}
                '''.format(gff=self.gff,
                           input=self.input().path,
                           output=self.output().path)


@inherits(GFFread)
class GetConsensusesWrapper(luigi.Task):
    lib_list = luigi.ListParameter()
    library = None
    consensus_type = None

    def requires(self):
        return {consensus_type: [self.clone_parent(library=library, consensus_type=consensus_type)
                                 for library in self.lib_list]
                for consensus_type in ['H1', 'H2', 'iupac-codes']}

    def output(self):
        return self.input()


@requires(GetConsensusesWrapper)
class GetAlignment(SlurmTask):
    min_cov = luigi.FloatParameter(default=0.8)
    min_indvs = luigi.FloatParameter(default=0.8)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return {'phy': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, self.output_prefix + ".phy")),}
                #'nex': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, self.output_prefix + ".nex"))}

    def work(self):
        import Bio
        import Bio.SeqIO
        import Bio.AlignIO
        import contextlib
        import numpy as np

        with contextlib.ExitStack() as stack, self.output()['phy'].open('w') as fphy:#, self.output()['nex'].open('w') as fnex:

            fhs = [stack.enter_context(open(fname.path)) for fname in self.input()['iupac-codes']]
            parsers = zip(*[Bio.SeqIO.parse(f, 'fasta') for f in fhs])
            msa = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''), id=lib) for lib in self.lib_list]

            for seqs in parsers:
                id, l = seqs[0].id, len(seqs[0])
                assert all([x.id == id for x in seqs]), "Fasta sequences not sorted!"

                coverage = 1 - np.array([x.seq.count('N') for x in seqs]) / l
                indvs = np.mean(coverage > self.min_cov)

                if indvs > self.min_indvs:
                    for (i, x) in enumerate(seqs):
                        # 3rd codon
                        msa[i] += x.seq[::3]

            Bio.AlignIO.write(Bio.Align.MultipleSeqAlignment(msa), fphy, 'phylip-relaxed')
            #Bio.AlignIO.write(Bio.Align.MultipleSeqAlignment(msa), fnex, 'nexus')


@requires(GetAlignment)
class RAxML(SlurmExecutableTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 10
        self.partition = "nbi-long,RG-Diane-Saunders"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mle', "RAxML_result." + self.output_prefix))

    def work_script(self):
        return '''#!/bin/bash
               source raxml-8.2.9;

               cd {output_dir}
               rm {output_dir}/RAxML*
               set -euo pipefail

               raxmlHPC-PTHREADS-SSE3 -T {n_cpu} -s {input} -m GTRGAMMA -n {suffix}.temp -p 100 ;

               mv RAxML_result.{suffix}.temp RAxML_result.{suffix}
               '''.format(output_dir=os.path.split(self.output().path)[0],
                          n_cpu=self.n_cpu,
                          input=self.input()['phy'].path,
                          suffix=self.output_prefix)


@requires(GetAlignment)
class RAxML_Bootstrap(SlurmExecutableTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 20
        self.partition = "nbi-long,RG-Diane-Saunders"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'bootstraps', "RAxML_bootstrap." + self.output_prefix))

    def work_script(self):
        return '''#!/bin/bash
               source raxml-8.2.9;
               cd {output_dir}
               rm {output_dir}/RAxML*
               set -euo pipefail

               raxmlHPC-PTHREADS-SSE3 -T {n_cpu} -s {input} -m GTRGAMMA -n {suffix}.temp -p 100 -b 1234 -N 100;

               mv RAxML_bootstrap.{suffix}.temp RAxML_bootstrap.{suffix}
               '''.format(output_dir=os.path.split(self.output().path)[0],
                          n_cpu=self.n_cpu,
                          input=self.input()['phy'].path,
                          suffix=self.output_prefix)


@inherits(RAxML)
@inherits(RAxML_Bootstrap)
class RAxML_Combine(SlurmExecutableTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 2
        self.partition = "nbi-medium"

    def requires(self):
        return {'mle': self.clone(RAxML),
                'bootstrap': self.clone(RAxML_Bootstrap)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'comb', "RAxML_bipartitionsBranchLabels." + self.output_prefix))

    def work_script(self):
        return '''#!/bin/bash
               source raxml-8.2.9;
               cd {output_dir}
               rm {output_dir}/RAxML*
               set -euo pipefail

               raxmlHPC-PTHREADS-SSE3 -T {n_cpu} -m GTRGAMMA -p 100 -f b -t {mle} -z {bootstrap} -n {suffix}.temp

               mv RAxML_bipartitionsBranchLabels.{suffix}.temp RAxML_bipartitionsBranchLabels.{suffix}
               '''.format(output_dir=os.path.split(self.output().path)[0],
                          n_cpu=self.n_cpu,
                          bootstrap=self.input()['bootstrap'].path,
                          mle=self.input()['mle'].path,
                          suffix=self.output_prefix)


if __name__ == '__main__':
    multiprocessing.set_start_method('forkserver')
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=os.path.basename(__file__))

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]

    luigi.run(['RAxML_Combine', '--output-prefix', name,
                                '--lib-list', json.dumps(lib_list),
                                '--star-genome', os.path.join(utils.reference_dir, 'genome'),
                                '--gff', os.path.join(utils.reference_dir, 'PST_genes_final.gff3'),
                                '--reference', os.path.join(utils.reference_dir, 'PST130_contigs.fasta'),
                                '--mask', os.path.join(utils.reference_dir, 'PST130_RNASeq_collapsed_exons.bed')] + sys.argv[2:])
