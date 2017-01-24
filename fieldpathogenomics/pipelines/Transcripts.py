import os
import sys
import json
import time

import luigi
from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile

from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from fieldpathogenomics.luigi.uv import UVExecutableTask
from fieldpathogenomics.utils import CheckTargetNonEmpty
from fieldpathogenomics.pipelines.Library import MarkDuplicates

import logging
logger = logging.getLogger('luigi-interface')
alloc_log = logging.getLogger('alloc_log')
alloc_log.setLevel(logging.DEBUG)


picard = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar"
python = "source /usr/users/ga004/buntingd/FP_dev/dev/bin/activate"

# Ugly hack
script_dir = os.path.join(os.path.split(
    os.path.split(__file__)[0])[0], 'scripts')
log_dir = os.path.join(os.path.split(
    os.path.split(os.path.split(__file__)[0])[0])[0], 'logs')
os.makedirs(log_dir, exist_ok=True)


# -----------------------------StringTie------------------------------- #


@requires(MarkDuplicates)
class StringTie(SlurmExecutableTask, CheckTargetNonEmpty):

    library = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(
        default="/tgac/scratch/buntingd/", significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, self.library, 'stringtie.gtf'))

    def work_script(self):
        return '''#!/bin/bash
        source stringtie-1.3.0;
        set -euo pipefail

        stringtie {input} -p {n_cpu} > {output}.temp

        mv {output}.temp {output}
        '''.format(input=self.input().path,
                   output=self.output().path,
                   n_cpu=self.n_cpu,
                   )


@inherits(StringTie)
class StringTieMerge(SlurmExecutableTask, CheckTargetNonEmpty):

    lib_list = luigi.ListParameter()
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def requires(self):
        return [self.clone(StringTie, library=lib) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'stringtie.gtf'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source stringtie-1.3.0;
        set -euo pipefail

        echo '{input}' > {temp}
        stringtie  -p {n_cpu} --merge {temp} > {output}.temp

        mv {output}.temp {output}
        '''.format(input="\n".join([x.path for x in self.input()]),
                   output=self.output().path,
                   temp=self.temp.path,
                   n_cpu=self.n_cpu,
                   )

# -----------------------------Cufflinks------------------------------- #


@requires(MarkDuplicates)
class Cufflinks(SlurmExecutableTask, CheckTargetNonEmpty):

    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(
        default="/tgac/scratch/buntingd/", significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, self.library, "cufflinks.gtf"))

    def work_script(self):
        return '''#!/bin/bash
        source cufflinks-2.2.1;
        set -euo pipefail

        cufflinks {input} --quiet -o {temp_dir}

        mv {temp_dir}/transcripts.gtf {output}
        '''.format(input=self.input().path,
                   output=self.output().path,
                   temp_dir=os.path.join(self.scratch_dir, self.library, 'cufflinks_temp'))


@inherits(Cufflinks)
class CuffMerge(SlurmExecutableTask, CheckTargetNonEmpty):

    lib_list = luigi.ListParameter()
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def requires(self):
        return [self.clone(Cufflinks, library=lib) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'cufflinks.gtf'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source cufflinks-2.2.1;
        source python-2.7.10;

        set -euo pipefail

        mkdir -p {out_dir}

        echo '{input}' > {temp}
        cd {out_dir}
        cuffmerge  -p {n_cpu} {temp} -o {out_dir}

        mv {out_dir}/merged_asm/merged.gtf {output}
        '''.format(input="\n".join([x.path for x in self.input()]),
                   output=self.output().path,
                   out_dir=os.path.join(self.scratch_dir, 'cufflinks_merge'),
                   temp=self.temp.path,
                   n_cpu=self.n_cpu,
                   )

# -----------------------------Trinity------------------------------- #


@inherits(MarkDuplicates)
class MergeBam(SlurmExecutableTask, CheckTargetNonEmpty):
    lib_list = luigi.ListParameter()
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 16000
        self.n_cpu = 3
        self.partition = "tgac-medium"

    def requires(self):
        return [self.clone(MarkDuplicates, library=lib) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, 'merged.bam'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source samtools-1.3;
        set -euo pipefail

        echo '{input}' > {temp}
        samtools merge -f  {output}.temp.bam -b {temp} --threads 2

        mv {output}.temp.bam {output}
        '''.format(input="\n".join([x.path for x in self.input()]),
                   output=self.output().path,
                   temp=self.temp.path)


@requires(MergeBam)
class Trinity(UVExecutableTask, CheckTargetNonEmpty):

    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 100000
        self.n_cpu = 6
        self.host = 'uv2k2'

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'Trinity-GG.fasta'))

    def work_script(self):
        return '''#!/bin/bash
        export PATH='/tgac/software/testing/bin/:/tgac/software/production/bin/:'$PATH

        source trinityrnaseq-2.3.2;
        set -euo pipefail

        mkdir -p /tgac/scratch/buntingd/trinity
        cd /tgac/scratch/buntingd/trinity

        Trinity --genome_guided_bam {input} --max_memory 100G \
                --genome_guided_max_intron 10000 \
                --output /tgac/scratch/buntingd/trinity \
                --CPU {n_cpu}
        mv /tgac/scratch/buntingd/trinity/Trinity-GG.fasta {output}
        '''.format(input=self.input().path,
                   output=self.output().path,
                   n_cpu=self.n_cpu)


@requires(Trinity)
class GMAP(CheckTargetNonEmpty, SlurmExecutableTask):
    gmap_reference_name = luigi.Parameter(default='PST130')
    gmap_reference_path = luigi.Parameter(default='/tgac/workarea/collaborators/saunderslab/FP_pipeline/reference/gmap')

    def __init__(self, *args, **kwargs):
        self.mem = 4000
        self.n_cpu = 1
        self.partition = 'tgac-medium'

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'trinity.gtf'))

    def work_script(self):
        return '''!#/bin/bash
                source gmap-20160923;
                set -euo pipefail

                gmap -n 0 -D {db_path} -d {db} {fasta} -f gff3_gene > {output}.temp
                mv {output}.temp {output}

                '''.format(db_path=self.gmap_reference_path,
                           df=self.gmap_reference_name,
                           output=self.outfile().path)

# -----------------------------Portcullis------------------------------- #


@requires(MergeBam)
class PortcullisPrep(SlurmExecutableTask):
    reference = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 16000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, 'portcullis_prep'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source portcullis-1.0.0_beta6;
        set -euo pipefail

        portcullis prep --output {output}_temp {reference} {input}

        mv {output}_temp {output}
        '''.format(input=self.input().path,
                   reference=self.reference,
                   output=self.output().path)


@requires(PortcullisPrep)
class PortcullisJunc(SlurmExecutableTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 24000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, 'portcullis_junc'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source portcullis-1.0.0_beta6;
        set -euo pipefail

        portcullis junc --output {output}_temp/portcullis \
                        --orientation FR \
                        --strandedness unstranded \
                        {input}

        mv {output}_temp {output}
        '''.format(input=self.input().path,
                   output=self.output().path)

# -----------------------------Strawberry------------------------------- #


@requires(MarkDuplicates)
class Strawberry(SlurmExecutableTask, CheckTargetNonEmpty):

    library = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(
        default="/tgac/scratch/buntingd/", significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, self.library, 'strawberry.gtf'))

    def work_script(self):
        return '''#!/bin/bash
        source strawberry-0.8.3;
        set -euo pipefail

        strawberry --output-dir {temp_dir} {input}

        mv {temp_dir}/assembled_transcripts.gtf {output}
        '''.format(input=self.input().path,
                   output=self.output().path,
                   temp_dir=os.path.join(self.scratch_dir, self.library, 'strawberry_temp'))


@requires(Strawberry)
class FixStrawberry(luigi.Task):

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, self.library, 'strawberry_fix.gtf'))

    def run(self):
        import re

        trans_regex = re.compile("^PST\S+")
        exon_regex = re.compile("^(?!P).+\tStrawberry\texon")

        with self.input().open() as fin:
            with self.output().open('w') as fout:
                contig = None
                for line in fin:
                    if line[0] == '#':
                        fout.write(line)
                        continue

                    m = trans_regex.match(line)
                    if m:
                        contig = line.split('\t')[0]
                        fout.write(line)
                    else:
                        fout.write(exon_regex.sub(contig + '\tStrawberry\texon', line))


@inherits(FixStrawberry)
class StrawberryMerge(SlurmExecutableTask, CheckTargetNonEmpty):

    lib_list = luigi.ListParameter()
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def requires(self):
        return [self.clone(FixStrawberry, library=lib) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'strawberry.gtf'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source cufflinks-2.2.1;
        source python-2.7.10;
        set -euo pipefail

        mkdir -p {out_dir}
        cd {out_dir}

        echo '{input}' > {temp}
        cuffmerge  -p {n_cpu} {temp} -o {out_dir}
        sed  -i 's/Cufflinks/strawberry/' {out_dir}/merged_asm/merged.gtf

        mv {out_dir}/merged_asm/merged.gtf {output}
        '''.format(input="\n".join([x.path for x in self.input()]),
                   output=self.output().path,
                   out_dir=os.path.join(self.scratch_dir, 'strawberry_merge'),
                   temp=self.temp.path,
                   n_cpu=self.n_cpu,
                   )

#   ----------------------------------------------------------------------#


@inherits(CuffMerge)
@inherits(StringTieMerge)
@inherits(PortcullisJunc)
@inherits(GMAP)
class TranscriptsWrapper(luigi.WrapperTask):

    def requires(self):
        yield self.clone(StringTieMerge)
        yield self.clone(CuffMerge)
        yield self.clone(Trinity)
        yield self.clone(PortcullisJunc)


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logging.disable(logging.DEBUG)
    timestr = time.strftime("%Y%m%d-%H%M%S")

    fh = logging.FileHandler(os.path.join(
        log_dir, os.path.basename(__file__) + "_" + timestr + ".log"))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    alloc_fh = logging.FileHandler(os.path.join(
        log_dir, os.path.basename(__file__) + "_" + timestr + ".salloc.log"))
    alloc_fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    alloc_fh.setFormatter(formatter)
    alloc_log.addHandler(alloc_fh)

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]

    luigi.run(['TranscriptsWrapper', '--lib-list', json.dumps(lib_list),
                                     '--star-genome', '/tgac/workarea/collaborators/saunderslab/Realignment/data/genome/',
                                     '--reference', '/tgac/workarea/collaborators/saunderslab/Realignment/data/PST130_contigs.fasta'] + sys.argv[2:])
