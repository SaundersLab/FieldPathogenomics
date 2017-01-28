import os
import sys
import json
import math

import luigi
from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile

from fieldpathogenomics.utils import CheckTargetNonEmpty
from fieldpathogenomics.SGUtils import ScatterBED, GatherVCF, ScatterVCF, GatherTSV
from fieldpathogenomics.luigi.scattergather import ScatterGather
import fieldpathogenomics.utils as utils
import fieldpathogenomics.picard.Library as Library

picard = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar"
gatk = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/gatk/3.6.0/x86_64/bin/GenomeAnalysisTK.jar "
snpeff = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/snpeff/4.3g/x86_64/snpEff.jar "
snpsift = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/snpeff/4.3g/x86_64/SnpSift.jar "

python = "source /usr/users/ga004/buntingd/FP_dev/dev/bin/activate"

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


@inherits(Library.HaplotypeCaller)
class gVCFs(luigi.Task, CheckTargetNonEmpty):
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)

    output_prefix = luigi.Parameter()
    reference = luigi.Parameter()
    lib_list = luigi.ListParameter()
    library = None

    def requires(self):
        return [self.clone(Library.HaplotypeCaller, library=lib) for lib in self.lib_list]

    def output(self):
        return self.input()


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
        return LocalTarget(os.path.join(self.scratch_dir, self.output_prefix, "combined", self.output_prefix + "_" + str(self.idx) + ".g.vcf"))

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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_raw.vcf.gz"))

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


@ScatterGather(ScatterVCF, GatherTSV, N_scatter)
@inherits(GenotypeGVCF)
class VariantsToTable(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def requires(self):
        return self.clone(GenotypeGVCF)

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, 'QC', self.output_prefix + ".tsv"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  set -eo pipefail

                  $gatk -T VariantsToTable -R {reference} --allowMissingData -F CHROM \
                                                                             -F POS \
                                                                             -F QD \
                                                                             -F FS \
                                                                             -F DP \
                                                                             -F QUAL \
                                                                             -GF GQ \
                                                                             -GF RGQ \
                                                                             -GF DP \
                                                                             -V {input} -o {output}.temp

                  mv {output}.temp {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem * self.n_cpu))


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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_filtered.vcf.gz"))

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
class GetSNPs(SlurmExecutableTask, CheckTargetNonEmpty):
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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_SNPs.vcf.gz"))

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


@inherits(GetSNPs)
class VCFtoHDF5(SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 16000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_SNPs.hd5"))

    def work_script(self):
        self.temp1 = TemporaryFile()
        cache_dir = os.path.join(self.scratch_dir, self.output_prefix, 'SNPs.hd5.cache')
        os.makedirs(cache_dir)
        return '''#!/bin/bash
                {python}
                set -eo pipefail
                gzip -cd {input} > {temp1}

                vcf2npy --vcf {temp1} --array-type calldata_2d --output-dir {cache_dir} --compress
                vcf2npy --vcf {temp1} --array-type variants --output-dir {cache_dir} --compress

                vcfnpy2hdf5 --vcf {temp1} --input-dir {cache_dir} --ouput {output}.temp

                mv {output}.temp {output}
                '''.format(python=python,
                           input=self.input().path,
                           temp1=self.temp1.path + '.vcf',
                           cache_dir=cache_dir,
                           output=self.output().path)


@requires(VcfToolsFilter)
class VariantsEval(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, 'QC', self.output_prefix + ".variant_eval"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  set -eo pipefail

                  $gatk -T VariantEval -R {reference} --eval {input} -o {output}.temp

                  mv {output}.temp {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem * self.n_cpu))


@requires(GetSNPs)
class SnpEff(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_SNPs_ann.vcf.gz"))

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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_SNPs_syn.vcf.gz"))

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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_INDELs_only.vcf.gz"))

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
class GetRefSNPs(SlurmExecutableTask, CheckTargetNonEmpty):
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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix + "_RefSNPs.vcf.gz"))

    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  source vcftools-0.1.13;
                  gatk='{gatk}'
                  set -eo pipefail

                  $gatk -T -T SelectVariants -V {input} -R {reference} \
                        --restrictAllelesTo BIALLELIC \
                        --selectTypeToInclude SYMBOLIC \
                        --selectTypeToInclude NO_VARIATION \
                        --selectTypeToInclude SNP \
                        --out {output}.temp.vcf.gz

                  # Filter out * which represents spanning deletions
                  gzip -cd {output}.temp.vcf.gz | grep -v $'\t\*\t' | bgzip -c > {output}.temp2.vcf.gz

                  mv {output}.temp2.vcf.gz {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem * self.n_cpu))

# ----------------------------------------------------------------------- #


@inherits(GetSyn)
@inherits(GetRefSNPs)
@inherits(GetINDELs)
@inherits(VariantsToTable)
@inherits(VariantsEval)
class CallsetWrapper(luigi.WrapperTask):

    def requires(self):
        yield self.clone(GetINDELs)
        yield self.clone(GetSyn)
        yield self.clone(VariantsToTable)
        yield self.clone(VariantsEval)
        yield self.clone(GetRefSNPs)


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=os.path.basename(__file__))

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]

    luigi.run(['CallsetWrapper', '--output-prefix', name,
                                 '--lib-list', json.dumps(lib_list),
                                 '--reference', '/tgac/workarea/collaborators/saunderslab/Realignment/data/PST130_contigs.fasta',
                                 '--mask', '/tgac/workarea/users/buntingd/realignment/PST130/Combined/PST130_RNASeq_collapsed_exons.bed'] + sys.argv[3:])
