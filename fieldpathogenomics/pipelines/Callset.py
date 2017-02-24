import os
import sys
import json
import math
import shutil

import luigi
from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile

import fieldpathogenomics
from fieldpathogenomics.utils import CheckTargetNonEmpty, gatk, snpeff, snpsift
from fieldpathogenomics.SGUtils import ScatterBED, GatherVCF, ScatterVCF
from fieldpathogenomics.luigi.scattergather import ScatterGather
from fieldpathogenomics.luigi.commit import CommittedTarget, CommittedTask
from fieldpathogenomics.luigi.notebook import NotebookTask
import fieldpathogenomics.utils as utils
import fieldpathogenomics.pipelines.Library as Library

python = "source /usr/users/ga004/buntingd/FP_dev/dev/bin/activate"

FILE_HASH = utils.file_hash(__file__)
VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]
PIPELINE = os.path.basename(__file__).split('.')[0]

'''
Guidelines for harmonious living:
--------------------------------
1. Tasks acting on fastq files should output() a list like [_R1.fastq, _R2.fastq]
2. Tasks acting on a bam should just output() a single Target
3. Tasks acting on a vcf should just output() a single Target

Notes
-------
job.mem is actually mem_per_cpu
'''

N_scatter = int(sys.argv[2]) if __name__ == '__main__' else 5


class GenomeContigs(luigi.ExternalTask):
    '''one per line list of contigs in the genome'''
    mask = luigi.Parameter()

    def output(self):
        return LocalTarget(self.mask)


@inherits(Library.CleanUpLib)
class gVCFs(luigi.Task, CheckTargetNonEmpty):
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)

    output_prefix = luigi.Parameter()
    reference = luigi.Parameter()
    lib_list = luigi.ListParameter()
    library = None

    def requires(self):
        return [self.clone(Library.CleanUpLib, library=lib) for lib in self.lib_list]

    def output(self):
        return [x['gvcf'] for x in self.input()]


@requires(gVCFs)
class CombineGVCFs(SlurmExecutableTask, CheckTargetNonEmpty):

    # Number of combined gVCFs to end up with
    N_gvcfs = luigi.IntParameter(default=5)
    idx = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 16000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.output_prefix, "combined", self.output_prefix + "_" + str(self.idx) + ".g.vcf"))

    def work_script(self):
        perfile = math.ceil(len(self.input()) / self.N_gvcfs)
        start_idx = perfile * self.idx
        end_idx = perfile * (self.idx + 1)
        self.variants = self.input()[
            start_idx:end_idx] if self.idx < self.N_gvcfs - 1 else self.input()[start_idx:]

        return '''#!/bin/bash
                source jre-8u92
                source gatk-3.6.0
                gatk='{gatk}'

                set -eo pipefail
                $gatk -T CombineGVCFs -R {reference} -o {output}.temp {variants}

                mv {output}.temp {output}
                '''.format(output=self.output().path,
                           gatk=gatk.format(mem=self.mem * self.n_cpu),
                           reference=self.reference,
                           variants="\\\n".join([" --variant " + lib.path for lib in self.variants]))


@inherits(CombineGVCFs)
class CombineGVCFsWrapper(luigi.Task):
    idx = None

    def requires(self):
        return [self.clone(CombineGVCFs, idx=idx) for idx in range(self.N_gvcfs)]

    def output(self):
        return self.input()


@ScatterGather(ScatterBED, GatherVCF, N_scatter)
@inherits(GenomeContigs)
@inherits(CombineGVCFsWrapper)
class GenotypeGVCF(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Combine the per sample g.vcfs into a complete callset
    :param str output_prefix: '''

    def requires(self):
        return [self.clone(GenomeContigs), self.clone(CombineGVCFsWrapper)]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, self.output_prefix + "_raw.vcf.gz"))

    def work_script(self):
        return '''#!/bin/bash
                source jre-8u92
                source gatk-3.6.0
                gatk='{gatk}'

                set -eo pipefail
                $gatk -T GenotypeGVCFs -R {reference} -L {intervals} -o {output}.temp.vcf.gz --includeNonVariantSites {variants}

                mv {output}.temp.vcf.gz {output}
                '''.format(output=self.output().path,
                           intervals=self.input()[0].path,
                           gatk=gatk.format(mem=self.mem * self.n_cpu),
                           reference=self.reference,
                           variants="\\\n".join([" --variant " + lib.path for lib in self.input()[1]]))


@ScatterGather(ScatterVCF, GatherVCF, N_scatter)
@inherits(GenotypeGVCF)
class VcfToolsFilter(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Applies hard filtering to the raw callset'''
    GQ = luigi.IntParameter(default=30)
    QD = luigi.IntParameter(default=5)
    FS = luigi.IntParameter(default=30)

    def requires(self):
        return self.clone(GenotypeGVCF)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 2
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, self.output_prefix + "_filtered.vcf.gz"))

    def work_script(self):
        self.temp1 = TemporaryFile()
        self.temp2 = TemporaryFile()

        return '''#!/bin/bash
                source vcftools-0.1.13;
                source bcftools-1.3.1;
                set -eo pipefail

                bcftools view --apply-filters . {input} -o {temp1} -O z --threads 1
                bcftools filter {temp1} -e "FMT/RGQ < {GQ} || FMT/GQ < {GQ} || QD < {QD} || FS > {FS}" --set-GTs . -o {temp2} -O z --threads 1
                vcftools --gzvcf {temp2} --recode --max-missing 0.000001 --stdout --bed {mask} | bgzip -c > {output}.temp

                mv {output}.temp {output}
                tabix -p vcf {output}
                '''.format(input=self.input().path,
                           output=self.output().path,
                           GQ=self.GQ,
                           QD=self.QD,
                           FS=self.FS,
                           mask=self.mask,
                           temp1=self.temp1.path,
                           temp2=self.temp2.path)


@ScatterGather(ScatterVCF, GatherVCF, N_scatter)
@inherits(VcfToolsFilter)
class GetSNPs(SlurmExecutableTask, CommittedTask, CheckTargetNonEmpty):
    '''Extracts just sites with only biallelic SNPs that have a least one variant isolate'''

    def requires(self):
        return self.clone(VcfToolsFilter)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return CommittedTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, self.output_prefix + "_SNPs.vcf.gz"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  source vcftools-0.1.13;
                  gatk='{gatk}'
                  set -eo pipefail

                  $gatk -T -T SelectVariants -V {input} -R {reference} --restrictAllelesTo BIALLELIC \
                                                                       --selectTypeToInclude SNP \
                                                                       --out {output}.temp.vcf.gz

                  # Filter out * which represents spanning deletions
                  gzip -cd {output}.temp.vcf.gz | grep -v $'\t\*\t' | bgzip -c > {output}.temp2.vcf.gz

                  rm {output}.temp.vcf.gz
                  mv {output}.temp2.vcf.gz {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem * self.n_cpu))


class VCFtoHDF5(SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 64000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(self.input().path.rsplit('.', maxsplit=1)[0] + ".hd5")

    def to_str_params(self, only_significant=False):
        sup = super().to_str_params(only_significant)
        extras = {'input': self.input().path}
        return dict(list(sup.items()) + list(extras.items()))

    def work_script(self):
        self.temp1 = TemporaryFile()
        cache_dir = self.output().path + '.cache'
        shutil.rmtree(cache_dir, ignore_errors=True)
        os.makedirs(cache_dir)
        return '''#!/bin/bash
                {python}
                rm {output}.temp
                set -eo pipefail
                gzip -cd {input} > {temp1}

                vcf2npy --vcf {temp1} --arity 'AD:6' --array-type calldata_2d --output-dir {cache_dir}
                vcf2npy --vcf {temp1} --arity 'AD:6' --array-type variants --output-dir {cache_dir}

                vcfnpy2hdf5 --vcf {temp1} --input-dir {cache_dir} --output {output}.temp

                mv {output}.temp {output}
                '''.format(python=python,
                           input=self.input().path,
                           temp1=self.temp1.path + '.vcf',
                           cache_dir=cache_dir,
                           output=self.output().path)


@inherits(GenotypeGVCF)
@inherits(VcfToolsFilter)
@inherits(GetSNPs)
class HD5s(luigi.WrapperTask):

    def requires(self):
        yield self.clone(requires(GenotypeGVCF)(VCFtoHDF5))
        yield self.clone(requires(VcfToolsFilter)(VCFtoHDF5))
        yield self.clone(requires(GetSNPs)(VCFtoHDF5))


@requires(requires(GetSNPs)(VCFtoHDF5))
@inherits(GetSNPs)
class SNPsNotebook(NotebookTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.notebook = os.path.join(utils.notebooks, 'Callset', 'SNPs.ipynb')
        self.vars_dict = {'SNPS_HD5': self.input().path}
        logger.info(str(self.vars_dict))

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, 'QC', 'SNPs.ipynb'))


@inherits(SNPsNotebook)
class QCNotebooks(luigi.WrapperTask):
    def requires(self):
        yield self.clone(SNPsNotebook)


@requires(GetSNPs)
class SnpEff(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, self.output_prefix + "_SNPs_ann.vcf.gz"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source vcftools-0.1.13;
                  snpeff='{snpeff}'
                  set -eo pipefail

                  $snpeff PST130 {input} | bgzip -c > {output}.temp

                  mv {output}.temp {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             snpeff=snpeff.format(mem=self.mem * self.n_cpu))


@requires(SnpEff)
class GetSyn(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, self.output_prefix + "_SNPs_syn.vcf.gz"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source vcftools-0.1.13;
                  snpsift='{snpsift}'
                  set -eo pipefail

                  $snpsift filter "ANN[*].EFFECT has 'synonymous_variant'" {input} | bgzip -c > {output}.temp

                  mv {output}.temp {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             snpsift=snpsift.format(mem=self.mem * self.n_cpu))


@ScatterGather(ScatterVCF, GatherVCF, N_scatter)
@inherits(VcfToolsFilter)
class GetINDELs(SlurmExecutableTask, CheckTargetNonEmpty):
    '''Get sites with MNPs'''

    def requires(self):
        return self.clone(VcfToolsFilter)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, self.output_prefix + "_INDELs_only.vcf.gz"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  set -eo pipefail

                  $gatk -T -T SelectVariants -V {input} -R {reference} --selectTypeToInclude MNP \
                                                                       --selectTypeToInclude MIXED \
                                                                       --out {output}.temp.vcf.gz

                  mv {output}.temp.vcf.gz {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem * self.n_cpu))


@ScatterGather(ScatterVCF, GatherVCF, N_scatter)
@inherits(VcfToolsFilter)
class GetRefSNPs(SlurmExecutableTask, CommittedTask, CheckTargetNonEmpty):
    '''Create a VCF with SNPs and include sites that are reference like in all samples'''

    def requires(self):
        return self.clone(VcfToolsFilter)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return CommittedTarget(os.path.join(self.base_dir, VERSION, PIPELINE, 'callsets', self.output_prefix, self.output_prefix + "_RefSNPs.vcf.gz"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  source vcftools-0.1.13;
                  gatk='{gatk}'
                  set -eo pipefail

                  $gatk -T SelectVariants -V {input} -R {reference} \
                        --selectTypeToInclude NO_VARIATION \
                        --selectTypeToInclude SNP \
                        --out {output}.temp.vcf.gz

                  # Filter out * which represents spanning deletions
                  gzip -cd {output}.temp.vcf.gz | grep -v $'[,\t]\*' | bgzip -c > {output}.temp2.vcf.gz

                  mv {output}.temp2.vcf.gz {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem * self.n_cpu))

# ----------------------------------------------------------------------- #


@inherits(GetSyn)
@inherits(GetRefSNPs)
@inherits(GetINDELs)
@inherits(HD5s)
@inherits(QCNotebooks)
class CallsetWrapper(luigi.WrapperTask):

    def requires(self):
        yield self.clone(GetINDELs)
        yield self.clone(GetSyn)
        yield self.clone(GetRefSNPs)
        yield self.clone(HD5s)
        yield self.clone(QCNotebooks)


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=os.path.basename(__file__))

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]

    luigi.run(['CallsetWrapper', '--output-prefix', name,
                                 '--lib-list', json.dumps(lib_list),
                                 '--star-genome', '/tgac/workarea/collaborators/saunderslab/Realignment/data/genome/',
                                 '--reference', '/tgac/workarea/collaborators/saunderslab/Realignment/data/PST130_contigs.fasta',
                                 '--mask', '/tgac/workarea/users/buntingd/realignment/PST130/Combined/PST130_RNASeq_collapsed_exons.bed'] + sys.argv[3:])
