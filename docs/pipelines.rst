Pipelines
===========

The entire Field Pathogenomics project is broken down into a series of Luigi pipelines. The Library pipeline is run by library and the Callset pipeline is run once per variant callset.

Library Pipeline
-----------------

The Library pipeline contains all the tasks that should be run on a per sample basis wrapped up in :class:`PerLibPipeline` and a batch task to run :class:`PerLibPipeline` on each library in lib-list

The key steps are:
    * Locating and concatenating FASTQ files
    * Read trimming for quality and TrueSeq adaptors using Trimommatic
    * Read QC with FastQC
    * Alignment to reference with STAR
    * Variant calling with HaplotypeCaller
    * Plotting allele frequency 
    
    
Inputs
^^^^^^^
**Reads** .fastq.gz

    The first task in the pipeline :class:`FetchFastqGZ` searches :param:`read-dir` for files matching the pattern `*library*_R1.fastq.gz` and `*library*_R2.fastq.gz`. 
    If multiple matches are found they are concatenated together in filename sorted order.
    Symbolic links are ignored
    
Outputs
^^^^^^^
**dedupped.bam**

    The aligned, cleaned, annotated de-duplicated BAM file and index.
    
**Log.final.out**

    The STAR alignment log 

**AlignmentStats**

    A row is inserted into the SQL table AlignmentStats containing the alignment statistics.
    
**LIBxxxxx.g.vcf**

    gVCF file and index output by HaplotypeCaller.
    
**QC/**

    Directory containing nucleotide distributions, quality scores and allele frequency distributions.
    
    
Flags
^^^^^

.. attribute:: --base-dir

    Base of the output directory structure:: 
        \|--BASE_DIR
            \|--libraries
                \|--LIBxxxxxx

.. attribute:: --scratch-dir

    Scratch directory

.. attribute:: --read-dir 

    File glob to pass to find to search for reads in. 

.. attribute:: --star-genome

    STAR genome index

.. attribute:: --snp-db 

    VCF of high quality a priori known SNPs to use for base quality score recalibration (BQSR).
    If none BQSR is not run

.. attribute:: --reference 

    Reference genome FASTA (must have index in same folder)

.. attribute:: --lib-list 

    JSON formatted list of librariries. The :py:class:`PerLibPipeline` is run once for each library in lib-list.


Typically the Library pipeline will be run in using a SLURM wrapper script similar to this

.. code-block:: bash

   #!/bin/sh
   #SBATCH --mem 20000 
   #SBATCH -n 1
   #SBATCH -N 1
    
    cd /tgac/workarea/collaborators/saunderslab/FP_pipeline
    source production/bin/activate
    export LUIGI_CONFIG_PATH=/tgac/workarea/collaborators/saunderslab/FP_pipeline/luigi.cfg

    srun python -m src.pipelines.Library all_libs.txt \
    --base-dir /tgac/workarea/collaborators/saunderslab/FP_pipeline/data/ \
    --scratch-dir /tgac/scratch/buntingd/ \
    --read-dir /tgac/data/read/*DianeSaunders* \
    --workers 250 

    

When running the :py:module:`Library` as  script like this it takes as its first argument a text file containing libraries one per line which is used to construct lib-list

Callset Pipeline
-----------------

The Callset pipeline take the individually called gVCF files produced by HaplotypeCaller and co-calls variants to produce the callset, which is then filtered and evaluated.

Key steps:
     * Use CombineGVCF to collapse the single sample gVCFs into a set of multisample gVCFs. Broad recommend ~200 samples per gVCF
     * Run GenotypeGVCF to call variants.
     * Filter with vcftools
     * Separate out SNP, INDELs and non-variant sites.
     
     TODO:
     * Evaluate quality of the callset and produce a report/write stats to SQL 
     
     
     
Inputs
^^^^^^^
**gVCF** .g.vcf 
    Individually called gVCFs output from Library pipeline
    For each library in --lib-list the pipeline will expect a gVCF --gvcf-dir/LIBxxxxxxx/LIBxxxxxxx.g.vcf

**Exon Mask** .bed 
    To speed up the analysis of RNAseq data use the existing gene annotations and only call variants in exon regions.
    
Outputs
^^^^^^^^
**_SNPs.vcf.gz**
    Filtered single nucleotide variants

**_RefSNPs.vcf.gz**
    Filtered single nucleotide variants and sites called as homokaryotic reference-like
**_INDELs.vcf.gz**
    Filtered multinucleotide variants
    
TODO:
    * QC plots 
    * Statistical summaries
    
    
    
Flags
^^^^^

.. attribute:: --output-prefix

    Name to attach to callset

.. attribute:: --mask
    BED file of regions on which to perform genotyping e.g. exons for RNAseq

.. attribute:: --base-dir

    Base of the output directory structure:: 
        \|--BASE_DIR
            \|--callsets
                \|--(output-prefix)

.. attribute:: --scratch-dir

    Scratch directory

.. attribute:: --gVCF-dir 

    Folder to look for gVCFs in
    
.. attribute:: --reference 

    Reference genome FASTA (must have index in same folder)

.. attribute:: --lib-list 

    JSON formatted list of librariries to jointly genotype


NB: To speed up the analysis this pipeline makes heavy use of ScatterGather.
It scatters the regions in --mask and passes these to GenotypeGVCF though the -L flag. It then scatters the resulting vcf.gz files line-by-line for filtering.
The number of scatter gather spits to make for this is set in Callset.py, if running as module as below this parameter is exposed as the second arg.

Like Library.py, Callset.py is designed to be run as a script taking a text file containing as list of libraries as the first arg and number if SG splits as the second.
The name of the libraries file is taken as the output-prefix for the callset.

 
.. code-block:: bash

    #!/bin/sh
    #SBATCH --mem 16000
    #SBTACH -N 1
    #SBATCH -n 1
    #SBATCH -c 1
    #SBATCH -p tgac-medium
    
    cd /tgac/workarea/collaborators/saunderslab/FP_pipeline/
    source /tgac/workarea/collaborators/saunderslab/FP_pipeline/production/bin/activate
    
    export LUIGI_CONFIG_PATH=/tgac/workarea/collaborators/saunderslab/FP_pipeline/luigi.cfg
    
    srun python -m src.pipelines.Callset all_libs.txt 50 \
    --base-dir /tgac/workarea/collaborators/saunderslab/FP_pipeline/data/ 
    --gVCF-dir /tgac/workarea/collaborators/saunderslab/FP_pipeline/data/libraries/
    --workers 150 

To run efficiently need 3x as many workers as scatter gather splits ie 50 and 150 here.

