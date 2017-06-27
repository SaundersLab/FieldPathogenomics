import os
import sys
import json

import luigi
from luigi import LocalTarget
from luigi.file import TemporaryFile

from bioluigi.slurm import SlurmExecutableTask, SlurmTask
from bioluigi.utils import CheckTargetNonEmpty
from bioluigi.decorators import requires, inherits

import fieldpathogenomics
import fieldpathogenomics.utils as utils
import fieldpathogenomics.pipelines.Library as Library

FILE_HASH = utils.file_hash(__file__)
PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = fieldpathogenomics.__version__.rsplit('.', 1)[0]

# -----------------------------StringTie------------------------------- #


@requires(Library.MarkDuplicates)
class StringTie(SlurmExecutableTask, CheckTargetNonEmpty):

    library = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, 'stringtie.gtf'))

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
    output_prefix = luigi.Parameter()
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def requires(self):
        return [self.clone(StringTie, library=lib) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'stringtie.gtf'))

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


@requires(Library.MarkDuplicates)
class Cufflinks(SlurmExecutableTask, CheckTargetNonEmpty):

    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(
        default="/tgac/scratch/buntingd/", significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "cufflinks.gtf"))

    def work_script(self):
        return '''#!/bin/bash
        source cufflinks-2.2.1;
        set -euo pipefail

        cufflinks {input} --quiet -o {temp_dir}

        mv {temp_dir}/transcripts.gtf {output}
        '''.format(input=self.input().path,
                   output=self.output().path,
                   temp_dir=os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, 'cufflinks_temp'))


@inherits(Cufflinks)
class CuffMerge(SlurmExecutableTask, CheckTargetNonEmpty):

    lib_list = luigi.ListParameter()
    library = None
    output_prefix = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def requires(self):
        return [self.clone(Cufflinks, library=lib) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'cufflinks_raw.gtf'))

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
                   out_dir=os.path.join(self.scratch_dir, VERSION, PIPELINE, 'cufflinks_merge'),
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
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'cufflinks.gtf'))

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


@inherits(Library.MarkDuplicates)
class MergeBam(SlurmExecutableTask, CheckTargetNonEmpty):
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)

    lib_list = luigi.ListParameter()
    output_prefix = luigi.Parameter()
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 16000
        self.n_cpu = 3
        self.partition = "nbi-medium"

    def requires(self):
        return [self.clone(Library.MarkDuplicates, library=lib) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.output_prefix, 'merged.bam'))

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
class Trinity(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 10
        self.partition = 'RG-Diane-Saunders,nbi-long'

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'Trinity-GG.fasta'))

    def work_script(self):
        return '''#!/bin/bash
        source trinityrnaseq-2.3.2;
        mkdir -p {scratch}
        set -euo pipefail

        cd {scratch}

        Trinity --genome_guided_bam {input} \
                --max_memory {mem} \
                --genome_guided_max_intron 10000 \
                --output {scratch} \
                --CPU {n_cpu}

        mv {scratch}/Trinity-GG.fasta {output}
        '''.format(input=self.input().path,
                   output=self.output().path,
                   n_cpu=self.n_cpu,
                   mem=int(0.95 * self.mem * self.n_cpu / 1000),
                   scratch=os.path.join(self.scratch_dir, VERSION, PIPELINE, self.output_prefix))


@requires(Trinity)
class GMAP(SlurmExecutableTask, CheckTargetNonEmpty):

    gmap_reference_name = luigi.Parameter(default='PST130')
    gmap_reference_path = luigi.Parameter(default='/tgac/workarea/collaborators/saunderslab/FP_pipeline/reference/gmap')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mem = 4000
        self.n_cpu = 1
        self.partition = 'nbi-medium'

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'trinity.gff3'))

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
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.output_prefix, 'portcullis_prep'))

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
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.output_prefix, 'portcullis_junc'))

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
        self.partition = "nbi-medium"

    def requires(self):
        return {'prep': self.clone(PortcullisPrep),
                'junc': self.clone(PortcullisJunc)}

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.output_prefix, 'portcullis_filter'))

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
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'configuration.yaml'))

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
                '''.format(python=utils.python,
                           list=self.list(),
                           temp=self.temp.path,
                           reference=self.reference,
                           portcullis=os.path.join(self.input()['portcullis'].path, 'portcullis.pass.junctions.bed'),
                           db=self.blast_db,
                           output=self.output().path)


@requires(MikadoConfigure)
class MikadoPrepare(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return {'gtf': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'mikado_prepared.gtf')),
                'fasta': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'mikado_prepared.fasta'))}

    def work_script(self):
        return '''#!/bin/bash
                  {python}
                  set -euo pipefail

                  mikado prepare  --strip_cds  -o {gtf}.temp -of {fasta}.temp --json-conf {conf} --log /dev/stderr

                  mv {gtf}.temp {gtf}
                  mv {fasta}.temp {fasta}
                '''.format(python=utils.python,
                           gtf=self.output()['gtf'].path,
                           fasta=self.output()['fasta'].path,
                           conf=self.input().path)


@requires(MikadoPrepare)
class BLAST(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 750
        self.n_cpu = 12
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'blast.xml.gz'))

    def work_script(self):
        return '''#!/bin/bash
                  source ncbi_blast+-2.2.30;
                  set -euo pipefail

                  blastx -max_target_seqs 5 \
                         -num_threads {n_cpu} \
                         -query {input} \
                         -outfmt 5 \
                         -db {blast_db} \
                         -evalue 0.000001 | sed '/^$/d' | gzip -c  > {output}.temp

                mv {output}.temp {output}
               '''.format(n_cpu=self.n_cpu,
                          input=self.input()['fasta'].path,
                          blast_db=self.blast_db,
                          output=self.output().path)


@requires(MikadoPrepare)
class TransDecoder(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'orfs.bed'))

    def work_script(self):
        return '''#!/bin/bash
                  source transdecoder-3.0.0
                  set -euo pipefail
                  mkdir -p {scratch}/transdecoder
                  cd {scratch}/transdecoder

                  TransDecoder.LongOrfs -t {input} 2> lo_log.txt
                  TransDecoder.Predict -t {input} --cpu {n_cpu}  2> pred_log.txt

                   mv {scratch}/transdecoder/{prefix}.transdecoder.bed {output}
               '''.format(input=self.input()['fasta'].path,
                          output=self.output().path,
                          prefix=os.path.split(self.input()['fasta'].path)[1],
                          scratch=os.path.join(self.scratch_dir, VERSION, PIPELINE),
                          n_cpu=self.n_cpu)


@inherits(MikadoConfigure)
@inherits(BLAST)
@inherits(TransDecoder)
class MikadoSerialise(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def requires(self):
        return {'conf': self.clone(MikadoConfigure),
                'blast': self.clone(BLAST),
                'orfs': self.clone(TransDecoder)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'mikado.db'))

    def work_script(self):
        return '''#!/bin/bash
                  {python}
                  rm {output}.temp
                  set -euo pipefail
                  cd {output_dir}

                  mikado serialise  --orfs {orfs} --xml {blast} --json-conf {conf} {output}.temp

                  mv {output}.temp {output}
                '''.format(python=utils.python,
                           output_dir=os.path.split(self.output().path)[0],
                           orfs=self.input()['orfs'].path,
                           blast=self.input()['blast'].path,
                           output=self.output().path,
                           conf=self.input()['conf'].path)


@inherits(MikadoConfigure)
@inherits(MikadoPrepare)
@inherits(MikadoSerialise)
class MikadoPick(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1500
        self.n_cpu = 8
        self.partition = "nbi-short"

    def requires(self):
        return {'conf': self.clone(MikadoConfigure),
                'db': self.clone(MikadoSerialise),
                'prep': self.clone(MikadoPrepare)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'mikado.gff'))

    def work_script(self):
        return '''#!/bin/bash
                  {python}
                  set -euo pipefail
                  cd {output_dir}

                  mikado pick -p {n_cpu} --json-conf {conf} -db {db} --loci_out {output}.temp --log /dev/stderr

                  mv {output}.temp.gff3 {output}
                '''.format(python=utils.python,
                           n_cpu=self.n_cpu,
                           output_dir=os.path.split(self.output().path)[0],
                           gff=self.input()['prep']['gtf'].path,
                           db=self.input()['db'].path,
                           output=self.output().path,
                           conf=self.input()['conf'].path)


@requires(MikadoPick)
class MikadoCompare(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1500
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'mikado', 'mikado.gff'))

    def work_script(self):
        return '''#!/bin/bash
                  {python}
                  set -euo pipefail
                  cd {output_dir}

                  mikado pick -p {n_cpu} --json-conf {conf} -db {db} --loci_out {output}.temp --log /dev/stderr

                  mv {output}.temp.gff3 {output}
                '''.format(python=utils.python,
                           n_cpu=self.n_cpu,
                           output_dir=os.path.split(self.output().path)[0],
                           gff=self.input()['prep']['gtf'].path,
                           db=self.input()['db'].path,
                           output=self.output().path,
                           conf=self.input()['conf'].path)

# ----------------------------------------------------------------------#


if __name__ == '__main__':
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=os.path.basename(__file__))

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]
    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]

    luigi.run(['MikadoPick', '--lib-list', json.dumps(lib_list),
                             '--output-prefix', name,
                             '--star-genome', os.path.join(utils.reference_dir, 'genome'),
                             '--reference', os.path.join(utils.reference_dir, 'PST130_contigs.fasta'),
                             '--blast-db', os.path.join(utils.reference_dir, '/uniprot/pst_uniprot.fasta')] + sys.argv[2:])
