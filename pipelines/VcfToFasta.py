## get_fasta

import os

import luigi
from luigi.contrib.slurm import SlurmExecutableTask
from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile

from fieldpathogenomics import GetRefSNPSs

class GetVCF(luigi.ExternalTask):
    '''Input VCF file, containing both SNPs and non-variant sites'''
    ref_snp_vcf = luigi.Parameter()
    base_dir = luigi.Parameter(default="/usr/users/ga004/buntingd/FP_dev/testing/", significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)

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
        return LocalTarget(os.path.join(self.scratch_dir, os.path.basename(self.ref_snp_vcf).split(".")[0]+".bcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                bcftools view {input} -o {output} -O b --threads 1
                '''.format(input=self.input().path,
                           output=self.output().path)


@requires(ConvertToBCF)
class GetSingleSample(SlurmExecutableTask):
    '''Pull a single sample out of the joint BCF file and use picard to compute a BED file where there is missing data'''
    lib = luigi.Parameter()
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 2
        self.partition = "tgac-short"

    def output(self):
        return[LocalTarget(os.path.join(self.scratch_dir, 'single_sample', self.lib + ".vcf.gz"),
               LocalTarget(os.path.join(self.scratch_dir, 'single_sample', self.lib + ".bed")]

    def work_script(self):
        
        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                source vcftools-0.1.13;
                picard='{picard}'
                
                bcftools view {input} -o {vcf} -O z -s {lib} --exclude-uncalled --no-update 
                tabix -p vcf {vcf}
                bgzip -dc {vcf} | $picard IntervalListTools I=/dev/stdin INVERT=true O=/dev/stdout | $picard IntervalListToBed I=/dev/stdin O={mask}
                
                '''.format(picard=picard.format(mem=self.mem),
                            input=self.input().path,
                            vcf=self.output()[0].path,
                            mask=self.output()[1],
                            lib=self.lib)


@requires(GetSingleSample)        
class BCFtoolsConsensus(SlurmExecutableTask):
    '''Apply the variants in the vcf file and mask with the missing data BED file.
       :param str consensus_type: can be H1 or H2 for applying one of the (pseudo) haplotypes or iupac-codes for ambiguous coding'''

    consensus_type = luigi.Parameter() # {'H1' , 'H2', 'iupac-codes'}
    reference = luigi.Parameter(default=/tgac/workarea/collaborators/saunderslab/Realignment/data/PST130_contigs.fasta)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-short"
    
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'fasta', 'genomes', self.lib+self.consensus_type+'.fasta' )
        
    def work_script(self):
        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                source vcftools-0.1.13;
                
                bcftools consensus {vcf} -f {reference} -s {lib} -m {mask} {consensus_type} > {output}
                samtools faidx {output}
                
                '''.format(vcf=self.output()[0].path,
                           mask=self.output()[1].path,
                           reference=self.reference,
                           lib=self.lib, 
                           output=self.output().path, 
                           consensus_type=self.consensus_type)

@requires(BCFtoolsConsensus)
class GFFread(SlurmExecutableTask):
    '''Pull out the spliced exons from the genomes'''
    gff = luig.Parameter()
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-short"
        
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'fasta', 'genes', self.lib+self.consensus_type+'.fasta' )
        
    def work_script(self):
        return '''#!/bin/bash -e
                source bcftools-1.3.1;
                source vcftools-0.1.13;
                
                gffread {gff} -g {input} -w /dev/stdout | fold -w 60 > {output} 
                
                '''.format(gff=self.gff,
                           input=self.input().path,
                           outputself.output().path)
                           
class CleanUp(SlurmExecutableTask):
    clean_dir = luigi.Parameter()
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 500
        self.n_cpu = 1
        self.partition = "tgac-short"
        
    def work_script(self):
        return '''#!/bin/bash -e
                rm -rf {clean_dir}'''.format(clean_dir=self.clean_dir)
                
    def complete(self):
        return not os.path.exists(self.clean_dir)


@inherits(GetVCF)
class GetConsensusesWrapper(luigi.WrapperTask):
    lib_list = luigi.ListParameter()        
    def requires(self):
        for lib in self.lib_list:
            for consensus_type in ['H1', 'H2', 'iupac-codes']:
                yield GFFread(lib=lib, consensus_type=consensus_type)
    
        yield CleanUp(clean_dir=os.path.join(self.scratch_dir, 'single_sample'))

                           
                           
                           
                           
