Using ScatterGather
===================

The scattergather decorator allows you to super easily make tasks parallel.
It breaks the input file into N_scatter separate files, runs the task on each on in parallel then joins the output together. For example GetSNPs goes through the vcf file and selects good snp sites. This task is trivially parallel because each line of the vcf file can be processed separately (except the header). The great thing is that this happens completely automagically, you can write the task exactly as you would for serial execution then add the decorator to make it parallel

.. code-block:: python

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
            self.partition = "nbi-medium"

        def output(self):
            return CommittedTarget(os.path.join(self.base_dir, VERSION, PIPELINE, self.output_prefix, self.output_prefix + "_SNPs.vcf.gz"))

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

The ScatterGather decorator requires three parameters: a scatter task, a gather task and a number to scatter.
The scatter task just describes how the input file is broken up and the gather task how the output file is reconstructed.
For example here ScatterVCF/GatherVCF needs to properly handle the vcf header rather than just leaving it on the first chunk.
There are a selection of Scatter/Gather task for various file types in SG_utils.

The only caveat is that the ScatterGather decorator cannot be used with the requires decorator, you must use inherit and manually define the requires method.


Additionally, if you have consecutive ScatterGather tasks the scheduler with not gather and the re-scatter the output of the first task, it will just pass straight to the input of the second task.


