import os
import sys
import json
import shutil
import sqlalchemy

import luigi
from luigi.contrib import sqla
from luigi.util import requires, inherits
from luigi import LocalTarget

import fieldpathogenomics
from fieldpathogenomics.luigi.slurm import SlurmExecutableTask, SlurmTask
from fieldpathogenomics.utils import CheckTargetNonEmpty, picard, gatk, trimmomatic
from fieldpathogenomics.luigi.commit import CommittedTarget, CommittedTask

import fieldpathogenomics.utils as utils

FILE_HASH = utils.file_hash(__file__)
PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]

'''

TODO: Migrate making the STAR reference to luigi and correctly set the genome column in AlginmentStats

Guidelines for harmonious living:
--------------------------------
1. Tasks acting on fastq files should output() a list like [_R1.fastq, _R2.fastq]
2. Tasks acting on a bam should just output() a single Target
3. Tasks acting on a vcf should just output() a single Target

Notes
-------
job.mem is actually mem_per_cpu
'''


class FetchFastqGZ(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Fetches and concatenate the fastq.gz files for ``library`` from the /reads/ server
     :param str library: library name  '''

    library = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)
    read_dir = luigi.Parameter(default="/tgac/data/reads/*DianeSaunders*", significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "raw_R1.fastq.gz")),
                LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "raw_R2.fastq.gz"))]

    def work_script(self):
        return '''#!/bin/bash -e
                  set -euo pipefail

                  find {read_dir} -name "*{library}*_R1.fastq.gz" -type f  -print | sort | xargs cat  > {R1}.temp
                  find {read_dir} -name "*{library}*_R2.fastq.gz" -type f  -print | sort | xargs cat  > {R2}.temp

                  mv {R1}.temp {R1}
                  mv {R2}.temp {R2}
                 '''.format(read_dir=self.read_dir,
                            library=self.library,
                            R1=self.output()[0].path,
                            R2=self.output()[1].path)


@requires(FetchFastqGZ)
class Trimmomatic(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "filtered_R1.fastq.gz")),
                LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "filtered_R2.fastq.gz")),
                LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, self.library + "_trimmomatic.txt"))]

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source trimmomatic-0.30
               set -euo pipefail

               cd {scratch_dir}
               trimmomatic='{trimmomatic}'
               $trimmomatic PE -threads 4 {R1_in} {R2_in} -baseout temp.fastq.gz \
               ILLUMINACLIP:{adapters}:2:30:10:4 SLIDINGWINDOW:4:20 MINLEN:50 \
               2>&1 | sed 's/raw_R1.fastq.gz/{library}.fastq.gz/' > {log}.temp

               mv temp_1P.fastq.gz {R1_out}
               mv temp_2P.fastq.gz {R2_out}
               mv {log}.temp {log}

                '''.format(scratch_dir=os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library),
                           trimmomatic=trimmomatic.format(
                               mem=self.mem * self.n_cpu),
                           log=self.output()[2].path,
                           library=self.library,
                           R1_in=self.input()[0].path,
                           R2_in=self.input()[1].path,
                           adapters='/tgac/software/testing/trimmomatic/0.30/x86_64/bin/adapters/TruSeq.cat.fa',
                           R1_out=self.output()[0].path,
                           R2_out=self.output()[1].path)


@requires(FetchFastqGZ)
class FastxQC(SlurmExecutableTask):
    '''Runs Fastx toolkit to plot the nucleotide and base call quality score distributions '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        working_dir = os.path.join(self.base_dir, VERSION, PIPELINE, self.library)
        return {'stats_R1': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R1_stats.txt")),
                'stats_R2': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R2_stats.txt")),
                'boxplot_R1': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R1_quality.png")),
                'boxplot_R2': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R2_quality.png")),
                'nt_dist_R1': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R1_nt_distr.png")),
                'nt_dist_R2': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R2_nt_distr.png")),
                }

    def work_script(self):
        return '''#!/bin/bash
        source fastx_toolkit-0.0.13.2
        set -euo pipefail

        gzip -cd {R1_in} | fastx_quality_stats -o {stats_R1} -Q33
        gzip -cd {R2_in} | fastx_quality_stats -o {stats_R2} -Q33

        fastq_quality_boxplot_graph.sh -i {stats_R1} -o {boxplot_R1}
        fastq_quality_boxplot_graph.sh -i {stats_R2} -o {boxplot_R2}

        fastx_nucleotide_distribution_graph.sh -i {stats_R1} -o {nt_dist_R1}
        fastx_nucleotide_distribution_graph.sh -i {stats_R2} -o {nt_dist_R2}

        '''.format(R1_in=self.input()[0].path,
                   R2_in=self.input()[1].path,
                   stats_R1=self.output()['stats_R1'].path,
                   stats_R2=self.output()['stats_R2'].path,
                   boxplot_R1=self.output()['boxplot_R1'].path,
                   boxplot_R2=self.output()['boxplot_R2'].path,
                   nt_dist_R1=self.output()['nt_dist_R1'].path,
                   nt_dist_R2=self.output()['nt_dist_R2'].path)


@requires(FetchFastqGZ)
class FastxTrimmer(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Uses FastxTrimmer to remove Illumina adaptors and barcodes'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "filtered_R1.fastq.gz")),
                LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "filtered_R2.fastq.gz"))]

    def work_script(self):
        return '''#!/bin/bash
        source fastx_toolkit-0.0.13.2
        set -euo pipefail

        gzip -cd {R1_in} | fastx_trimmer -f14 -Q33 | gzip > {R1_out}.temp ;
        gzip -cd {R2_in} | fastx_trimmer -f14 -Q33 | gzip > {R2_out}.temp ;

        mv {R1_out}.temp {R1_out}
        mv {R2_out}.temp {R2_out}
        '''.format(R1_in=self.input()[0].path,
                   R2_in=self.input()[1].path,
                   R1_out=self.output()[0].path,
                   R2_out=self.output()[1].path)


@requires(Trimmomatic)
class Star(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Runs STAR to align to the reference :param str star_genome:'''
    star_genome = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 3000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def output(self):
        return {
            'star_bam': LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, 'Aligned.sortedByCoord.out.bam')),
            'star_log': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, 'Log.final.out'))
        }

    def work_script(self):
        return '''#!/bin/bash
                  source star-2.5.0a
                  set -euo pipefail

                  mkdir -p {scratch_dir}/star_temp
                  cd  {scratch_dir}/star_temp

                  STAR  --genomeDir {star_genome} \
                        --outSAMstrandField intronMotif \
                        --outSAMtype BAM SortedByCoordinate \
                        --runThreadN {n_cpu} \
                        --readFilesCommand gunzip -c \
                        --readFilesIn {R1} {R2}

                  mv {scratch_dir}/star_temp/Log.final.out {star_log}
                  mv {scratch_dir}/star_temp/Aligned.sortedByCoord.out.bam {star_bam}

                  '''.format(star_bam=self.output()['star_bam'].path,
                             star_log=self.output()['star_log'].path,
                             scratch_dir=os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library),
                             star_genome=self.star_genome,
                             n_cpu=self.n_cpu,
                             R1=self.input()[0].path,
                             R2=self.input()[1].path,)


@requires(Trimmomatic)
class FastQC(CheckTargetNonEmpty, SlurmExecutableTask):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Set the SLURM request params for this task
            self.mem = 2000
            self.n_cpu = 1
            self.partition = "tgac-medium"

        def output(self):
            return [LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, 'QC', 'R1', 'fastqc_data.txt')),
                    LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, 'QC', 'R2', 'fastqc_data.txt'))]

        def work_script(self):
            return '''#!/bin/bash
                    source fastqc-0.11.4
                    mkdir -p {output_dir}
                    set -euo pipefail

                    fastqc {R1_in} {R2_in} -o {output_dir} -t 1

                    cd {output_dir}
                    unzip filtered_R1_fastqc.zip
                    sed 's/Filename\tfiltered_R1.fastq.gz/Filename\t{lib}_R1/'  filtered_R1_fastqc/fastqc_data.txt > {R1_out}.temp

                    unzip filtered_R2_fastqc.zip
                    sed 's/Filename\tfiltered_R2.fastq.gz/Filename\t{lib}_R2/'  filtered_R2_fastqc/fastqc_data.txt > {R2_out}.temp

                    mv {R1_out}.temp {R1_out}
                    mv {R2_out}.temp {R2_out}
                    '''.format(output_dir=os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, 'FastQC'),
                               R1_in=self.input()[0].path,
                               R2_in=self.input()[1].path,
                               lib=self.library,
                               R1_out=self.output()[0].path,
                               R2_out=self.output()[1].path)


@requires(Star)
class AlignmentStats(sqla.CopyToTable):
    columns = [
        (["Library", sqlalchemy.String(64)], {}),
        (["input_reads", sqlalchemy.INTEGER], {}),
        (["input_len", sqlalchemy.FLOAT], {}),
        (["mapped_reads", sqlalchemy.INTEGER], {}),
        (["mapped_reads_pc", sqlalchemy.String(10)], {}),
        (["mapped_len", sqlalchemy.FLOAT], {}),
        (["mismatch_pc", sqlalchemy.String(10)], {}),
        (["datetime", sqlalchemy.String(25)], {}),
        (["genome", sqlalchemy.String(25)], {}),
        (["git_commit", sqlalchemy.String(40)], {}),
        (["pipeline_hash", sqlalchemy.String(40)], {}),
    ]

    star_keys = {
        "Library": "Library",
        "input_reads": 'Number of input reads',
        "input_len": 'Average input read length',
        "mapped_reads": 'Uniquely mapped reads number',
        "mapped_reads_pc": 'Uniquely mapped reads %',
        "mapped_len": 'Average mapped length',
        "mismatch_pc": 'Mismatch rate per base, %',
        "datetime": "Started job on"
    }

    connection_string = "mysql+pymysql://tgac:tgac_bioinf@tgac-db1.hpccluster/buntingd_fieldpathogenomics"
    table = "AlignmentStats"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        git_commit = utils.current_commit_hash(os.path.split(__file__)[0])
        pipeline_hash = utils.hash_pipeline(self)
        genome = os.path.split(os.path.dirname(self.star_genome))[1]
        star_log = utils.parseStarLog(
            self.input()['star_log'].path, self.library)

        self._rows = [[star_log[AlignmentStats.star_keys[x[0][0]]] for x in AlignmentStats.columns[
            :len(star_log)]] + [genome, git_commit, pipeline_hash]]

    def rows(self):
        return self._rows

    def update_id(self):
        return hash(str(self._rows))


@requires(Star)
class CleanSam(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1500
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, 'Aligned.out_cleaned.bam'))

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source picardtools-2.1.1
               picard='{picard}'
               set -euo pipefail

               $picard CleanSam VERBOSITY=ERROR QUIET=true I={input} O={output}.temp

               mv {output}.temp {output}
                '''.format(input=self.input()['star_bam'].path,
                           output=self.output().path,
                           picard=picard.format(mem=self.mem * self.n_cpu))


@requires(CleanSam)
class AddReadGroups(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Sets the read group to the sample name, required for GATK'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 3000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, 'rg_added_sorted.bam'))

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source picardtools-2.1.1
               picard='{picard}'
               set -euo pipefail

               $picard AddOrReplaceReadGroups VERBOSITY=ERROR QUIET=true I={input} O={output}.temp SO=coordinate RGID=Star RGLB={lib} RGPL=Ilumina RGPU=Ilumina RGSM={lib}

               mv {output}.temp {output}
                '''.format(input=self.input().path,
                           output=self.output().path,
                           lib=self.library,
                           picard=picard.format(mem=self.mem * self.n_cpu))


@requires(AddReadGroups)
class MarkDuplicates(CheckTargetNonEmpty, CommittedTask, SlurmExecutableTask):
    '''Marks optical/PCR duplicates'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return CommittedTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, self.library + '.bam'))

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source picardtools-2.1.1
               picard='{picard}'
               set -euo pipefail

               $picard MarkDuplicates VERBOSITY=ERROR QUIET=true I={input} O={output}.temp CREATE_INDEX=false VALIDATION_STRINGENCY=SILENT M=/dev/null

               mv {output}.temp {output}
                '''.format(input=self.input().path,
                           output=self.output().path,
                           picard=picard.format(mem=self.mem * self.n_cpu))


@requires(MarkDuplicates)
class PortcullisFilterBam(SlurmExecutableTask):
    '''Removes reads correspoding to incorrect splice junctions'''

    portcullis_junc = luigi.Parameter(default='')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        if self.portcullis_junc == '':
            return self.input()
        else:
            return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, 'portcullis.bam'))

    def on_success(self):
        if self.portcullis_junc == '':
            luigi.Task.on_success(self)
        else:
            SlurmExecutableTask.on_success(self)

    def on_failure(self, e):
        if self.portcullis_junc == '':
            luigi.Task.on_failure(self, e)
        else:
            SlurmExecutableTask.on_failure(self, e)

    def run(self):
        if self.portcullis_junc != '':
            super().run()

    def work_script(self):
        return '''#!/bin/bash
               source portcullis-1.0.0_beta6;
               set -euo pipefail

               portcullis bamfilt  --output {output}.temp {junc} {bam}
               mv {output}.temp {output}
                '''.format(junc=self.portcullis_junc,
                           bam=self.input().path,
                           output=self.output().path)


@requires(PortcullisFilterBam)
class BaseQualityScoreRecalibration(SlurmExecutableTask):
    '''Runs BQSR. Because this requires a set of high quality SNPs to use
    as a ground truth we bootstrap this by first running the pipeline without
    BQSR then running again using the best SNPs of the first run.

    This is achieved by conditionally overriding run() on whether a snp_db is given
    '''
    snp_db = luigi.Parameter(default='')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        if self.snp_db == '':
            return self.input()
        else:
            return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, 'recalibrated.bam'))

    def on_success(self):
        if self.snp_db == '':
            luigi.Task.on_success(self)
        else:
            SlurmExecutableTask.on_success(self)

    def on_failure(self, e):
        if self.snp_db == '':
            luigi.Task.on_failure(self, e)
        else:
            SlurmExecutableTask.on_failure(self, e)

    def run(self):
        if self.snp_db != '':
            super().run()

    def work_script(self):
        recal = os.path.join(self.base_dir, VERSION, PIPELINE, self.library, self.library + "_recal.tsv")
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  set -euo pipefail

                  $gatk -T BaseRecalibrator  -R {reference}  -I {input}  -knownSites {snp_db}  -o {recal}
                  $gatk -T PrintReads -R {reference} -I {input} -BQSR {recal} -o {output}.temp

                  mv {output}.temp {output}
                '''.format(gatk=gatk.format(mem=self.mem * self.n_cpu),
                           input=self.input().path,
                           output=self.output().path,
                           reference=self.reference,
                           recal=recal)


@requires(BaseQualityScoreRecalibration)
class SplitNCigarReads(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Required by GATK, breaks up reads spanning introns'''
    reference = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, 'split.bam'))

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source gatk-3.6.0
               gatk='{gatk}'
               picard='{picard}'
               set -euo pipefail

               $picard BuildBamIndex VERBOSITY=ERROR QUIET=true I={input}
               $gatk -T SplitNCigarReads --logging_level ERROR -R {reference} -I {input} -o {output}.temp -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

               mv {output}.temp.bai {output}.bai
               mv {output}.temp {output}
                '''.format(input=self.input().path,
                           output=self.output().path,
                           picard=picard.format(mem=self.mem * self.n_cpu),
                           gatk=gatk.format(mem=self.mem * self.n_cpu),
                           reference=self.reference)


@requires(SplitNCigarReads)
class HaplotypeCaller(CheckTargetNonEmpty, CommittedTask, SlurmExecutableTask):
    '''Per sample SNP calling'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 1
        self.partition = "tgac-medium"
        self.sbatch_args = '--constraint=intel'


    def output(self):
        return CommittedTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.library, self.library + ".g.vcf"))

    def work_script(self):
        return '''#!/bin/bash
                source jre-8u92
                gatk='{gatk}'
                set -euo pipefail

                $gatk -T HaplotypeCaller  \
                      -R {reference} \
                      -I {input} \
                      -dontUseSoftClippedBases\
                      --emitRefConfidence GVCF \
                      -o {output}.temp.g.vcf

                mv {output}.temp.g.vcf {output}
        '''.format(input=self.input().path,
                   output=self.output().path,
                   gatk=gatk.format(mem=self.mem * self.n_cpu),
                   reference=self.reference)


@requires(HaplotypeCaller)
class PlotAlleleFreq(SlurmTask):
    '''Make plots of the ranked allele frequencies to identify mixed isolates'''

    DP_thresh = luigi.IntParameter(default=10)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, self.library, 'QC', self.library + "_allele_freqs.pdf"))

    def work(self):
        import vcfnp
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        variants = vcfnp.variants(self.input().path)
        calldata_2d = vcfnp.calldata_2d(self.input().path)

        var = np.logical_and(variants['ALT'] != b'<NON_REF>', variants['DP'] > self.DP_thresh)
        counts = np.sort(calldata_2d['AD'][var][:, 0, :], axis=1)[:, ::-1]
        freqs = counts / calldata_2d['DP'][var]
        third = 1 - np.sum(freqs, axis=1)

        df = pd.DataFrame(np.hstack((freqs, third.reshape((-1, 1)))))
        df[df == 0] = float('nan')

        df.hist(sharex=True, sharey=True, range=(0, 1), bins=20)
        plt.gcf().suptitle(self.library, fontsize=20)
        plt.gcf().text(0.5, 0.04, 'Allele frequency', ha='center')
        plt.gcf().text(0.02, 0.5, 'Counts', va='center', rotation='vertical')
        plt.gcf().savefig(self.output().path)


@inherits(Trimmomatic)
@inherits(FastxQC)
@inherits(FastQC)
class CombinedQC(luigi.WrapperTask):
    '''Wrapper task that runs all the QC type tasks library'''

    def requires(self):
        yield self.clone(FastQC)
        yield self.clone(FastxQC)


@inherits(MarkDuplicates)
@inherits(CombinedQC)
@inherits(HaplotypeCaller)
class PerLibPipeline(luigi.WrapperTask):
    '''Wrapper task that runs all tasks on a single library'''

    def requires(self):
        return {'qc': self.clone(CombinedQC),
                'gvcf': self.clone(HaplotypeCaller),
                'bam': self.clone(MarkDuplicates)}

    def output(self):
        return self.input()


@requires(PerLibPipeline)
class CleanUpLib(luigi.Task):
    priority = 100

    def run(self):
        shutil.rmtree(os.path.join(self.scratch_dir, VERSION, PIPELINE,
                                   self.library), ignore_errors=True)

    def complete(self):
        return self.clone_parent().complete() and not os.path.exists(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library))

    def output(self):
        return self.input()


@inherits(CleanUpLib)
class LibraryBatchWrapper(luigi.WrapperTask):
    '''Wrapper task to execute the per library part of the pipline on all
        libraries in :param list lib_list:'''
    lib_list = luigi.ListParameter()
    # This is a bit of a hack, it allows us to pass parameters to LibraryBatchWrapper and have them propagate
    # down to all calls to PerLibPipeline.
    library = None

    def requires(self):
        return [self.clone_parent(library=lib.rstrip()) for lib in self.lib_list]

    def output(self):
        return self.input()


# ----------------------------------------------------------------------- #


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=PIPELINE)

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    luigi.run(['LibraryBatchWrapper',
               '--lib-list', json.dumps(lib_list),
               '--star-genome', '/nbi/Research-Groups/JIC/Diane-Saunders/FP_pipeline/reference/genome/',
               '--reference', '/nbi/Research-Groups/JIC/Diane-Saunders/FP_pipeline/reference/PST130_contigs.fasta'] + sys.argv[2:])
