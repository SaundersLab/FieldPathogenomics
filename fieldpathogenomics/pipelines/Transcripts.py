import os
import sys
import json
import time

import luigi
from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile

from fieldpathogenomics.luigi.slurm import SlurmExecutableTask, SlurmTask
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
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'cufflinks_raw.gtf'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source cufflinks-2.2.1;
        source python-2.7.10;

        set -euo pipefail

        mkdir -p {out_dir}

        echo '{input}' > {temp}
        cd {out_dir}
        cuffmerge  -p {n_cpu}  -o {out_dir}  {temp}

        mv {out_dir}/merged.gtf {output}
        '''.format(input="\n".join([x.path for x in self.input()]),
                   output=self.output().path,
                   out_dir=os.path.join(self.scratch_dir, 'cufflinks_merge'),
                   temp=self.temp.path,
                   n_cpu=self.n_cpu,)


@requires(CuffMerge)
class AddTranscripts(SlurmTask):
    '''The gtf file produced by cuffmerge has no transcript features, not sure why??!
        This task reconstructs the transcript features use the transcript_id tag of the exons
        and taking the start/stop of the first/last exons with a given transcript_id is the start/stop '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'cufflinks.gtf'))

    def work(self):
        import re
        with self.input().open() as fin, self.output().open('w') as fout:
            current_tid = None
            tid = re.compile('transcript_id "(\S+)";')
            gid = re.compile('gene_id "(\S+)";')
            exons, start, stop = [], float('inf'), 0

            for line in fin:
                exon_start, exon_stop = line.split('\t')[3:5]
                tran = tid.search(line).groups()[0]

                if tran != current_tid:
                    if current_tid is not None:

                        # Flush previous transcript
                        template = exons[0].split('\t')
                        gene = gid.search(exons[0]).groups()[0]
                        fout.write('\t'.join(template[:2]) +
                                   '\ttranscript\t{0}\t{1}\t.\t{2}\t.\tgene_id "{3}"; transcript_id "{4}";\n'.format(
                                   start, stop, template[6], gene, current_tid))
                        fout.writelines(exons)

                    # And start new one
                    current_tid = tran
                    exons, start, stop = [], float('inf'), 0

                exons.append(line)
                start = min(start, int(exon_start))
                stop = max(stop, int(exon_stop))

            # Final flush
            template = exons[0].split('\t')
            gene = gid.search(exons[0]).groups()[0]
            fout.write('\t'.join(template[:2]) +
                       '\ttranscript\t{0}\t{1}\t.\t{2}\t.\tgene_id "{3}"; transcript_id "{4}";\n'.format(
                       start, stop, template[6], gene, current_tid))
            fout.writelines(exons)


# -----------------------------Trinity------------------------------- #


@inherits(MarkDuplicates)
class MergeBam(SlurmExecutableTask, CheckTargetNonEmpty):
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)

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
class GMAP(SlurmExecutableTask, CheckTargetNonEmpty):

    gmap_reference_name = luigi.Parameter(default='PST130')
    gmap_reference_path = luigi.Parameter(default='/tgac/workarea/collaborators/saunderslab/FP_pipeline/reference/gmap')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mem = 4000
        self.n_cpu = 1
        self.partition = 'tgac-medium'

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'trinity.gff3'))

    def work_script(self):
        return '''#!/bin/bash
                source gmap-20160923;
                set -euo pipefail

                gmap -n 0 -D {db_path} -d {db} {fasta} -f gff3_gene > {output}.temp
                mv {output}.temp {output}

                '''.format(db_path=self.gmap_reference_path,
                           db=self.gmap_reference_name,
                           output=self.output().path,
                           fasta=self.input().path)

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


@inherits(PortcullisJunc)
@inherits(PortcullisPrep)
class PortcullisFilter(SlurmExecutableTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def requires(self):
        return {'prep': self.clone(PortcullisPrep),
                'junc': self.clone(PortcullisJunc)}

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, 'portcullis_filter'))

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
        source portcullis-1.0.0_beta6;
        set -euo pipefail

        portcullis filter --output {output}_temp/portcullis {prep} {tab}

        mv {output}_temp {output}
        '''.format(prep=self.input()['prep'].path,
                   tab=os.path.join(self.input()['junc'].path, 'portcullis.junctions.tab'),
                   output=self.output().path)

# ----------------------------------------------------------------------#


@inherits(AddTranscripts)
@inherits(StringTieMerge)
@inherits(PortcullisFilter)
@inherits(GMAP)
class TranscriptsWrapper(luigi.Task):

    def requires(self):
        return {'stringtie': self.clone(StringTieMerge),
                'cufflinks': self.clone(AddTranscripts),
                'trinity': self.clone(GMAP),
                'portcullis': self.clone(PortcullisFilter)}

    def output(self):
        return self.input()

# -----------------------------Mikado------------------------------- #


@requires(TranscriptsWrapper)
class MikadoConfigure(SlurmExecutableTask, CheckTargetNonEmpty):

    blast_db = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'transcripts', 'mikado', 'configuration.yaml'))

    def list(self):
        '''Mikado requires a table of the transciptome assemblies to use.
           Construct this from the task input'''
        l = []
        for k, v in self.input().items():
            if k == 'cufflinks':
                l.append(v.path + "\tcu\tTrue")
            elif k == 'stringtie':
                l.append(v.path + "\tst\tTrue")
            elif k == 'class':
                l.append(v.path + "\tcl\tTrue")
            elif k == 'trinity':
                l.append(v.path + "\ttr\tFalse")
        return "\n".join(l)

    def work_script(self):
        self.temp = TemporaryFile()
        return '''#!/bin/bash
                   {python}
                   set -euo pipefail

                   echo '{list}' > {temp}

                   mikado configure --list {temp} \
                                    --reference {reference} \
                                    --mode permissive \
                                    --scoring plants.yaml  \
                                    --junctions {portcullis} \
                                    -bt {db} \
                                    {output}.temp
                mv {output}.temp {output}
                '''.format(python=python,
                           list=self.list(),
                           temp=self.temp.path,
                           reference=self.reference,
                           portcullis=os.path.join(self.input()['portcullis'].path, 'portcullis.pass.junctions.bed'),
                           db=self.blast_db,
                           output=self.output().path)


# ----------------------------------------------------------------------#

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
                                     '--reference', '/tgac/workarea/collaborators/saunderslab/Realignment/data/PST130_contigs.fasta',
                                     '--blast_db', '/tgac/workarea/collaborators/saunderslab/FP_pipeline/reference/uniprot.fasta'] + sys.argv[2:])
