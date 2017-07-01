Writing a new task
===================

Slurm Executable Task
---------------------

The bioluigi.SlurmExecutableTask is the most common task type used in the fieldpathogenomics pipelines, it just runs an external executable in a SLURM job. As an example I'll use running Fastx Trimmer, the complete task looks like this:

    .. code-block:: python

        @requires(FetchFastqGZ)
        class FastxTrimmer(CheckTargetNonEmpty, SlurmExecutableTask):
            '''Uses FastxTrimmer to remove Illumina adaptors and barcodes'''

            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)

                self.mem = 1000
                self.n_cpu = 1
                self.partition = "nbi-medium"

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


The first line
    .. code-block:: python

        @requires(c)

sets up the input requirements for the task. The output() method of the task given here become the input() method of this task. In this case FetchFastqGZ().output() returns a list containing the raw R1 and R2 fastq.gz files. Using the requires decorator also handles propagating luigi parameters down the task graph. For example FetchFastqGZ has a parameter library and so this task also inherits the library parameter.


The next line
    .. code-block:: python

        class FastxTrimmer(CheckTargetNonEmpty, SlurmExecutableTask):

defines the task name, which must be unique within the pipeline. It also sets the task type as SlurmExecutableTask and add the mixin CheckTargetNonEmpty, which checks that the output file produced is not just an empty file.

The next section

    .. code-block:: python

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            self.mem = 1000
            self.n_cpu = 1
            self.partition = "nbi-medium"

defines the slurm parameters for the task. **NB: mem is actually mem-per-cpu**. It's also crucial to make sure the super().__init__ line is there.

The output method

    .. code-block:: python

        def output(self):
            return [LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "filtered_R1.fastq.gz")),
                    LocalTarget(os.path.join(self.scratch_dir, VERSION, PIPELINE, self.library, "filtered_R2.fastq.gz"))]


tells the scheduler what files to expect the task to produce. The actual string path of a LocalTarget is accessed through LocalTarget().path.
Note the 'VERSION, PIPELINE' structure, this is important to make sure different versions of the same pipeline can't overwrite each others' data.

Finally them actual script for the task is defined in work_script, this method should return the bash script as a string.

    .. code-block:: python

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

The source command is required by the cluster to set up the program paths. Using 'set -euo pipefail' is important to make sure the task will fail if any part of the bash script fails, however the cluster source script has some (harmless) bash errors in that will cause the task to fail so has to come *before* set -euo!

In the command itself regular string substitution is used to put in the paths of the input and output file, the same can be done for parameters etc.

The final ' mv {R1_out}.temp {R1_out}' is very important for luigi. Luigi decides whether a task is complete or not by checking to see whether its output files exist or not. The problem is if a task dies half way though leaving a half written output file, in this case we want it to be marked as incomplete and run again but because the output exists luigi thinks it's done. The solution is to write the output to a temporary file and then at the very last second move it to the final destination.

Slurm Task
------------

The bioluigi.SlurmTask class allows you to write python code inline in the pipeline but have it execute in a separate slurm job.

    .. code-block:: python

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
                        'nex': LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, self.output_prefix + ".nex"))}

            def work(self):
                import Bio
                import Bio.SeqIO
                import Bio.AlignIO
                import contextlib
                import numpy as np

                with contextlib.ExitStack() as stack, self.output()['phy'].open('w') as fphy, self.output()['nex'].open('w') as fnex:

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
                    Bio.AlignIO.write(Bio.Align.MultipleSeqAlignment(msa), fnex, 'nexus')

The main structure is virtually identical to that of SlurmExectuableTask except instead of writing a work_script() method you write a work() method.
When the task is executed the whole task object is pickled and transported to the slurm work node where the work() then called. This means you have access to the task object as self from the work code.
Note that the general interpreter scope is not copied, so we have to begin by importing modules.
Because work() has access to self it can use self.input() and self.output() to access files. In this case no final mv is necessary because by using the self.output().open('w') construct (rather than open(self.output().path, 'w')) the final atomic move is performed automagically for you by luigi when the context is closed.

Notebook Tasks
--------------

I really like Jupyter notebooks as a way of working with data! The notebook task takes a template notebook and specialises it by populating special variables then runs it as part of the pipeline. For explain notebooks are used to QC callsets, there is a standard template of the analyses to be run that is specialised with the path of the input vcf.

.. code-block:: python

    @requires(HD5s)
    class SNPsNotebook(NotebookTask):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.mem = 8000
            self.n_cpu = 1
            self.partition = "nbi-medium"

            self.notebook = os.path.join(utils.notebooks, 'Callset', 'SNPs.ipynb')
            self.vars_dict = {'SNPS_HD5': self.input()['snps'].path}
            logger.info(str(self.vars_dict))

        def output(self):
            return LocalTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, 'QC', 'SNPs.ipynb'))


The self.notebook property points to the location of the template notebook, which is stored under version control.
The vars_dict property defines the variables to be substituted. In the template notebook there is a cell that looks like

.. code-block:: python

    ##luigi-vars
    SNPS_HD5 = 'default'

When the task is run this is replace by the values of SNPS_HD5 in vars_dict.
The completed notebook is output to output() as both a .ipynb copy (for editing) and a html copy (for viewing)

Database Tasks
--------------

Using the luigi module it is easy to write pipeline data to an SQL database. This is standard and is explained in the luigi docs. A real example of this is from the genome assembly pipeline, this task runs abyss-fac to calculate assembly N50 and store it in a table

.. code-block:: python

    class AbyssFac(sqla.CopyToTable):
        columns = [
            (["Task", sqlalchemy.String(20)], {}),
            (["K", sqlalchemy.INTEGER], {}),
            (["soap_k", sqlalchemy.INTEGER], {}),
            (["n", sqlalchemy.FLOAT], {}),
            (["n:500", sqlalchemy.FLOAT], {}),
            (["L50", sqlalchemy.FLOAT], {}),
            (["min", sqlalchemy.FLOAT], {}),
            (["N80", sqlalchemy.FLOAT], {}),
            (["N50", sqlalchemy.FLOAT], {}),
            (["N20", sqlalchemy.FLOAT], {}),
            (["E-size", sqlalchemy.FLOAT], {}),
            (["max", sqlalchemy.FLOAT], {}),
            (["sum", sqlalchemy.FLOAT], {}),
            (["path", sqlalchemy.String(500)], {})
        ]

        connection_string = "mysql+pymysql://tgac:tgac_bioinf@tgac-db1.hpccluster/buntingd_pstgenome"
        table = "abyssfac"

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def get_abyss(self):
            r = subprocess.run("source abyss-1.9.0; abyss-fac " + self.input().path,
                               stdout=subprocess.PIPE, shell=True, universal_newlines=True)
            r.check_returncode()
            values = r.stdout.split("\n")[1].split('\t')
            return [float(x) for x in values[:-1]] + [values[-1]]

        def rows(self):
            abyss = self.get_abyss()
            try:
                soap_k = self.soap_k
            except AttributeError:
                soap_k = -1
            self._rows = [[self.get_task_family()[:20]] + [self.K, soap_k] + abyss]
            return self._rows


The important thing to note here is the connection_string as this is specialised to the TGAC cluster. Also, this is run by the worker process so don't do anything to computationally intensive! Make a separate task to do calculations, output them to a file then make the CopyToTable task just read the file!

Committed Tasks
---------------

Committed tasks were something I experimented with but didn't really use. The idea is for important files used in downstream analyses like alignments and snp calls we need to be able to verify that the haven't been corrupted or edited since they were created. To do this the CommittedTask computes a checksum of the CommittedTarget after the task completes and stores this in a database along with the git commit hash of the code

.. code-block:: python

    @requires(SplitNCigarReads)
    class HaplotypeCaller(CheckTargetNonEmpty, CommittedTask, SlurmExecutableTask):
        '''Per sample SNP calling'''

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Set the SLURM request params for this task
            self.mem = 6000
            self.n_cpu = 1
            self.partition = "nbi-long,RG-Diane-Saunders"
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

