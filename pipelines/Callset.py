import os,sys, json,shutil
import multiprocessing
from time import sleep
import time
import logging
logger = logging.getLogger('luigi-interface')
alloc_log = logging.getLogger('alloc_log')
alloc_log.setLevel(logging.DEBUG)

import luigi
from luigi.contrib.slurm import SlurmExecutableTask
from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile
from src.utils import CheckTargetNonEmpty

from Library import LibraryBatchWrapper

picard="java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar"
gatk="java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/gatk/3.6.0/x86_64/bin/GenomeAnalysisTK.jar "
trimmomatic="java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/trimmomatic/0.36/x86_64/bin/trimmomatic-0.36.jar "
python="source /usr/users/ga004/buntingd/FP_dev/dev/bin/activate"

# Ugly hack
script_dir = os.path.join(os.path.split(os.path.split(__file__)[0])[0], 'scripts')
log_dir = os.path.join(os.path.split(os.path.split(os.path.split(__file__)[0])[0])[0], 'logs')
os.makedirs(log_dir, exist_ok=True)

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

@requires(LibraryBatchWrapper)        
class GenotypeGVCF(SlurmExecutableTask):
    '''Combine the per sample g.vcfs into a complete callset
    :param str output_prefix: '''
    output_prefix = luigi.Parameter(default="genotypes")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 32000
        self.n_cpu = 1
        self.partition = "tgac-medium"
    
    def input(self):
        ## THIS IS A MASSIVE FAT HACK
        return  [LocalTarget(os.path.join(self.base_dir, 'libraries', library, library + ".g.vcf")) for library in self.lib_list]
        
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix+"_raw.vcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash -e
                source jre-8u92
                source gatk-3.6.0
                gatk='{gatk}'
                $gatk -T GenotypeGVCFs -R {reference} -o {temp_out} --includeNonVariantSites {variants}
                
                
                mv {temp_out} {output}
                '''.format(temp_out=os.path.join(os.path.split(self.output().path)[0],"temp.vcf.gz"),
                           output=self.output().path,
                           gatk=gatk.format(mem=self.mem*self.n_cpu),
                           reference=self.reference,
                           variants="\n".join(["--variant "+ lib.path +" \\" for lib in self.input()]) )

@requires(GenotypeGVCF)
class ScatterVCF(SlurmExecutableTask):
    N_scatter = luigi.IntParameter(default=5)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"
    
    def run(self):
        if self.N_scatter > 1:
            super().run()
        elif self.N_scatter < 1:
            raise Exception("N_scatter must be >0")
    
    def output(self):
        self.dir = os.path.join(self.scratch_dir, self.output_prefix)
        if self.N_scatter > 1:
            return [LocalTarget(os.path.join(self.dir,'raw_{0}.vcf.gz'.format(i))) for i in range(self.N_scatter)]
        elif self.N_scatter == 1:
            return [self.input()]
        
    def work_script(self):
        return '''#!/bin/bash -e
                source vcftools-0.1.13;
                mkdir -p {dir}/temp
                
                bgzip -cd {input} | python {script_dir}/spilt_VCF.py {dir}/temp/raw {N_scatter}
                
                mv {dir}/temp/* {dir}
                rmdir {dir}/temp
                '''.format(dir=self.dir,
                           script_dir=script_dir,
                           input=self.input().path,
                           N_scatter=self.N_scatter)

@requires(ScatterVCF)
class VcfToolsFilter(SlurmExecutableTask):
    '''Applies hard filtering to the raw callset'''
    GQ = luigi.IntParameter(default=30)
    QD = luigi.IntParameter(default=5)
    FS = luigi.IntParameter(default=30)
    mask = luigi.Parameter(default="PST130_RNASeq_collapsed_exons.bed")
    scatter_idx = luigi.IntParameter()
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 2
        self.partition = "tgac-medium"
        
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_filtered"+str(self.scatter_idx)+".vcf.gz"))
    
    def work_script(self):
        self.temp1 = TemporaryFile()
        self.temp2 = TemporaryFile()
        
        return '''#!/bin/bash -e
                source vcftools-0.1.13;
                source bcftools-1.3.1;
                
                bcftools view --apply-filters . {input} -o {temp1} -O z --threads 1
                bcftools filter {temp1} -e "FMT/RGQ < {GQ} || FMT/GQ < {GQ} || QD < {QD} || FS > {FS}" --set-GTs . -o {temp2} -O z --threads 1
                vcftools --gzvcf {temp2} --recode --max-missing 0.000001 --stdout --bed {mask} | bgzip -c > {output}.temp
                
                mv {output}.temp {output}
                tabix -p vcf {output}
                '''.format(input=self.input()[self.scatter_idx].path,
                           output=self.output().path,
                           GQ=self.GQ,
                           QD=self.QD,
                           FS=self.FS,
                           mask=self.mask,
                           temp1=self.temp1.path,
                           temp2=self.temp2.path)

@requires(VcfToolsFilter)
class GetSNPs(SlurmExecutableTask):
    '''Extracts just sites with only biallelic SNPs that have a least one variant isolate'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"
    
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_SNPs_only"+str(self.scatter_idx)+".vcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash -e
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  $gatk -T -T SelectVariants -V {input} -R {reference} --restrictAllelesTo BIALLELIC --selectTypeToInclude SNP --out /dev/stdout | gzip -c > {output}.temp
                  
                  mv {output}.temp {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem*self.n_cpu))

@requires(VcfToolsFilter)
class GetINDELs(SlurmExecutableTask):
    '''Get sites with MNPs'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"
    
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_INDELs_only"+str(self.scatter_idx)+".vcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash -e
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  $gatk -T -T SelectVariants -V {input} -R {reference} --selectTypeToInclude MNP  --selectTypeToInclude MIXED   --out /dev/stdout | gzip -c > {output}.temp
                  
                  mv {output}.temp {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem*self.n_cpu))

@requires(VcfToolsFilter)
class GetRefSNPSs(SlurmExecutableTask):
    '''Create a VCF with SNPs and include sites that are reference like in all samples'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-medium"
    
    def output(self):
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_RefSNPs"+str(self.scatter_idx)+".vcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash -e
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  $gatk -T -T SelectVariants -V {input} -R {reference} --restrictAllelesTo BIALLELIC --selectTypeToInclude SYMBOLIC --selectTypeToInclude NO_VARIATION  --selectTypeToInclude SNP  --out /dev/stdout | \
                  gzip -c > {output}.temp
                  
                  mv {output}.temp {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem*self.n_cpu))

@inherits(GetSNPs)
class GatherSNPs(SlurmExecutableTask):
    scatter_idx=None
    def requires(self):
        return [self.clone_parent(scatter_idx=i) for i in range(self.N_scatter)]
    def output(self):
        return  LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_SNPs.vcf.gz"))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"
        
    def work_script(self):
        return '''#!/bin/bash -e
                picard='{picard}'
                source vcftools-0.1.13;
                source jre-8u92
                
                $picard MergeVcfs O={output}.temp {in_flags} 
                
                mv {output}.temp.vcf.gz {output}
                '''.format(picard=picard.format(mem=self.mem*self.n_cpu),
                           output=self.output().path,
                           in_flags="I=" + " I=".join([x.path for x in self.input()])
                           )

@inherits(GetRefSNPSs)
class GatherRefSNPs(SlurmExecutableTask):
    scatter_idx=None
    def requires(self):
        return [self.clone_parent(scatter_idx=i) for i in range(self.N_scatter)]
    def output(self):
        return  LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_RefSNPs.vcf.gz"))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"
        
    def work_script(self):
        return '''#!/bin/bash -e
                picard='{picard}'
                source vcftools-0.1.13;
                source jre-8u92
                
                $picard MergeVcfs O={output}.temp {in_flags} 
                
                mv {output}.temp.vcf.gz {output}
                '''.format(picard=picard.format(mem=self.mem*self.n_cpu),
                           output=self.output().path,
                           in_flags="I=" + " I=".join([x.path for x in self.input()])
                           )

@inherits(GetINDELs)
class GatherINDELs(SlurmExecutableTask):
    scatter_idx=None
    def requires(self):
        return [self.clone_parent(scatter_idx=i) for i in range(self.N_scatter)]
    def output(self):
        return  LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_INDELs.vcf.gz"))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"
        
    def work_script(self):
        return '''#!/bin/bash -e
                picard='{picard}'
                source vcftools-0.1.13;
                source jre-8u92
                
                $picard MergeVcfs O={output}.temp.vcf.gz {in_flags} 
                
                mv {output}.temp.vcf.gz {output}
                '''.format(picard=picard.format(mem=self.mem*self.n_cpu),
                           output=self.output().path,
                           in_flags="I=" + " I=".join([x.path for x in self.input()])
                           )

#-----------------------------------------------------------------------#

@inherits(GatherSNPs)
@inherits(GatherRefSNPs)
@inherits(GatherINDELs)
class SnpCalling(luigi.WrapperTask):
    def requires(self):
        yield self.clone(GatherINDELs)
        yield self.clone(GatherSNPs)
        yield self.clone(GatherRefSNPs)


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logging.disable(logging.DEBUG)
    timestr = time.strftime("%Y%m%d-%H%M%S")
    
    fh = logging.FileHandler(os.path.join(log_dir, os.path.basename(__file__) + "_" + timestr + ".log"))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    alloc_fh = logging.FileHandler(os.path.join(log_dir, os.path.basename(__file__) + "_" + timestr + ".salloc.log"))
    alloc_fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    alloc_fh.setFormatter(formatter)
    alloc_log.addHandler(alloc_fh)
    
    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file]
    
    luigi.run(['SnpCalling', '--lib-list', json.dumps(lib_list),
                             '--star-genome', '/tgac/workarea/collaborators/saunderslab/Realignment/data/genome/',
                             '--reference', '/tgac/workarea/collaborators/saunderslab/Realignment/data/PST130_contigs.fasta',
                             '--mask', '/tgac/workarea/users/buntingd/realignment/PST130/Combined/PST130_RNASeq_collapsed_exons.bed'] + sys.argv[2:])
