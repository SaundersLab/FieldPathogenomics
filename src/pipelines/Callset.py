import os,sys, json
import time,math
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
from src.SGUtils import ScatterBED, GatherVCF, ScatterVCF
from src.luigi.scattergather import ScatterGather

picard="java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar"
gatk="java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/gatk/3.6.0/x86_64/bin/GenomeAnalysisTK.jar "
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

N_scatter = int(sys.argv[2]) if __name__ == '__main__' else 5

class GenomeContigs(luigi.ExternalTask):
    '''one per line list of contigs in the genome'''
    mask = luigi.Parameter()
    def output(self):
        return LocalTarget(self.mask)

class gVCFs(luigi.ExternalTask):
    '''For each lib in lib_list get its gVCF'''
    lib_list = luigi.ListParameter()
    gVCF_dir = luigi.Parameter()
    
    output_prefix = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)
    reference = luigi.Parameter()
    
    def output(self):
        return [LocalTarget(os.path.join(self.gVCF_dir, lib, lib+'.g.vcf')) for lib in self.lib_list]

@requires(gVCFs)
class CombineGVCFs(SlurmExecutableTask, CheckTargetNonEmpty):
    
    N_gvcfs = luigi.IntParameter(default=5) # Number of combined gVCFs to end up with
    idx = luigi.IntParameter()
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 16000
        self.n_cpu = 1
        self.partition = "tgac-medium"
        
    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, self.output_prefix, "combined", self.output_prefix +"_" + str(self.idx) + ".g.vcf"))
    
    def work_script(self):
        perfile = math.ceil(len(self.input())/self.N_gvcfs)    
        start_idx = perfile*self.idx
        end_idx = perfile*(self.idx+1)
        self.variants = self.input()[start_idx:end_idx] if self.idx < self.N_gvcfs-1 else self.input()[start_idx:]
        
        return '''#!/bin/bash
                source jre-8u92
                source gatk-3.6.0
                gatk='{gatk}'
                
                set -eo pipefail
                $gatk -T CombineGVCFs -R {reference} -o {output}.temp {variants}
                
                mv {output}.temp {output}
                '''.format(output=self.output().path,
                           gatk=gatk.format(mem=self.mem*self.n_cpu),
                           reference=self.reference,
                           variants="\\\n".join([" --variant "+ lib.path for lib in self.variants]) )

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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix, self.output_prefix+"_raw.vcf.gz"))
    
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
                           gatk=gatk.format(mem=self.mem*self.n_cpu),
                           reference=self.reference,
                           variants="\\\n".join([" --variant "+ lib.path for lib in self.input()[1:]]) )

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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_filtered.vcf.gz"))
    
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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_SNPs.vcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  set -eo pipefail
                  
                  $gatk -T -T SelectVariants -V {input} -R {reference} --restrictAllelesTo BIALLELIC --selectTypeToInclude SNP --out {output}.temp.vcf.gz
                  
                  mv {output}.temp.vcf.gz {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem*self.n_cpu))

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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_INDELs_only.vcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  set -eo pipefail
                  
                  $gatk -T -T SelectVariants -V {input} -R {reference} --selectTypeToInclude MNP  --selectTypeToInclude MIXED   --out {output}.temp.vcf.gz
                  
                  mv {output}.temp.vcf.gz {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem*self.n_cpu))

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
        return LocalTarget(os.path.join(self.base_dir, 'callsets', self.output_prefix , self.output_prefix + "_RefSNPs.vcf.gz"))
        
    def work_script(self):
        return '''#!/bin/bash
                  source jre-8u92
                  source gatk-3.6.0
                  gatk='{gatk}'
                  set -eo pipefail
                  
                  $gatk -T -T SelectVariants -V {input} -R {reference} --restrictAllelesTo BIALLELIC --selectTypeToInclude SYMBOLIC --selectTypeToInclude NO_VARIATION  --selectTypeToInclude SNP  --out {output}.temp.vcf.gz
                  
                  mv {output}.temp.vcf.gz {output}
                  '''.format(input=self.input().path,
                             output=self.output().path,
                             reference=self.reference,
                             gatk=gatk.format(mem=self.mem*self.n_cpu))

#-----------------------------------------------------------------------#

@inherits(GetSNPs)
@inherits(GetRefSNPs)
@inherits(GetINDELs)
class CallsetWrapper(luigi.WrapperTask):
    def requires(self):
        yield self.clone(GetINDELs)
        yield self.clone(GetSNPs)
        yield self.clone(GetRefSNPs)

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

    name = os.path.split(sys.argv[1])[1].split('.', 1)[0]
    
    luigi.run(['CallsetWrapper', '--output-prefix', name,
                                 '--lib-list', json.dumps(lib_list),
                                 '--reference', '/tgac/workarea/collaborators/saunderslab/Realignment/data/PST130_contigs.fasta',
                                 '--mask', '/tgac/workarea/users/buntingd/realignment/PST130/Combined/PST130_RNASeq_collapsed_exons.bed'] + sys.argv[3:])
