# get_fasta

import os
import sys
import json
import shutil
import logging

import luigi
from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile

from SNP_calling import GetRefSNPSs, picard


class GetVCF(luigi.ExternalTask):
    '''Input VCF file, containing both SNPs and non-variant sites'''
    ref_snp_vcf = luigi.Parameter()
    base_dir = luigi.Parameter(
        default="/usr/users/ga004/buntingd/FP_dev/testing/", significant=False)
    scratch_dir = luigi.Parameter(
        default="/tgac/scratch/buntingd/", significant=False)

    def output(self):
        return LocalTarget(self.ref_snp_vcf)


@requires(GetVCF)
class ConvertToBCF(SlurmExecutableTask):
    '''Use bcftools view to convert the vcf to bcf, its worth doing this conversion 
     as the bcf formatted file is much faster for separating the samples than the vcf'''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 2
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, os.path.basename(self.ref_snp_vcf).split(".")[0] + ".bcf.gz"))

    def work_script(self):
        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                bcftools view {input} -o {output} -O b --threads 1
                '''.format(input=self.input().path,
                           output=self.output().path)


@requires(ConvertToBCF)
class GetSingleSample(SlurmExecutableTask):
    '''Pull a single sample out of the joint BCF file and use picard to compute a BED file where there is missing data'''
    library = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 2
        self.partition = "tgac-short"

    def output(self):
        return[LocalTarget(os.path.join(self.scratch_dir, 'single_sample', self.library + ".vcf.gz")),
               LocalTarget(os.path.join(self.scratch_dir, 'single_sample', self.library + ".bed"))]

    def work_script(self):

        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                source vcftools-0.1.13;
               source jre-8u92
               source picardtools-2.1.1
               picard='{picard}'
                
                bcftools view {input} -o {vcf} -O z -s {library} --exclude-uncalled --no-update 
                tabix -p vcf {vcf}
                $picard IntervalListTools I={vcf} INVERT=true O=/dev/stdout | $picard IntervalListToBed I=/dev/stdin O={mask}
                
                '''.format(picard=picard.format(mem=self.mem),
                           input=self.input().path,
                           vcf=self.output()[0].path,
                           mask=self.output()[1].path,
                           library=self.library)


@requires(GetSingleSample)
class BCFtoolsConsensus(SlurmExecutableTask):
    '''Apply the variants in the vcf file and mask with the missing data BED file.
       :param str consensus_type: can be H1 or H2 for applying one of the (pseudo) haplotypes or iupac-codes for ambiguous coding'''

    consensus_type = luigi.Parameter()  # {'H1' , 'H2', 'iupac-codes'}
    reference = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'fasta', 'genomes', self.library + "_" + self.consensus_type + '.fasta'))

    def work_script(self):
        ct_flag = "--" + self.consensus_type if self.consensus_type == 'iupac-codes' else "-" + \
            self.consensus_type
        return '''#!/bin/bash -e
        
                source bcftools-1.3.1;                
                bcftools consensus {vcf} -f {reference} -s {library} -m {mask} {consensus_type} > {output}
                
                source samtools-1.3;
                samtools faidx {output}
                
                '''.format(vcf=self.input()[0].path,
                           mask=self.input()[1].path,
                           reference=self.reference,
                           library=self.library,
                           output=self.output().path,
                           consensus_type=ct_flag)


@requires(BCFtoolsConsensus)
class GFFread(SlurmExecutableTask):
    '''Pull out the spliced exons from the genomes'''
    gff = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'fasta', 'genes', self.library + self.consensus_type + '.fasta'))

    def work_script(self):
        return '''#!/bin/bash -e
                source gffread-0.9.8;
                
                gffread {gff} -g {input} -w /dev/stdout | fold -w 60 > {output} 
                
                '''.format(gff=self.gff,
                           input=self.input().path,
                           output=self.output().path)


@inherits(GFFread)
class GetConsensusesWrapper(luigi.WrapperTask):
    lib_list = luigi.ListParameter()

    def requires(self):
        for library in self.lib_list:
            for consensus_type in ['H1', 'H2', 'iupac-codes']:
                yield self.clone_parent(library=library, consensus_type=consensus_type)

    def on_success(self):
        '''If the task successfully completes clean up the temporary files'''
        shutil.rmtree(os.path.join(self.scratch_dir, 'single_sample'))
        return super().on_success()
# Same hack as in LibraryBatchWrapper, allows us to propagate all
# arguments of GFFread upwards apart from library
GetConsensusesWrapper.library = None
GetConsensusesWrapper.consensus_type = None

if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logging.disable(logging.DEBUG)
    with open(sys.argv[1], 'r') as librarys_file:
        lib_list = [line.rstrip() for line in librarys_file]

    luigi.run(['GetConsensusesWrapper', '--lib-list', json.dumps(lib_list),
                                        '--gff', '/usr/users/ga004/buntingd/FP_dev/testing/data/PST_genes_final.gff3',
                                        '--ref-snp-vcf', '/usr/users/ga004/buntingd/FP_dev/testing/callsets/genotypes/genotypes_RefSNPs.vcf.gz',
                                        '--reference', '/usr/users/ga004/buntingd/FP_dev/testing/data/PST130_contigs.fasta',
                                        '--workers', '3', ])
