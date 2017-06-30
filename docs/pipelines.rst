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

    The first task in the pipeline :class:`FetchFastqGZ` searches `read-dir` for files matching the pattern `*library*_R1.fastq.gz` and `*library*_R2.fastq.gz`.
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


Typically the Library pipeline will be run indirect as part of the Callset pipeline, though could be run directly to just do QC/alignments.

When running the :py:module:`Library` as  script like this it takes as its first argument a text file containing libraries one per line which is used to construct lib-list

Callset Pipeline
-----------------

The Callset pipeline take the individually called gVCF files produced by HaplotypeCaller and co-calls variants to produce the callset, which is then filtered and evaluated.

Key steps:
     * Use CombineGVCF to collapse the single sample gVCFs into a set of multisample gVCFs. Broad recommend ~200 samples per gVCF
     * Run GenotypeGVCF to call variants.
     * Filter with vcftools
     * Separate out SNP, INDELs and non-variant sites.
     * Use SNPeff to get synonymous variant
     * Store the callset in HD5 format for easy access
     * Compute QC statistics and output to Jupyter notebook reports


Inputs
^^^^^^^
**gVCF** .g.vcf
    Individually called gVCFs output from Library pipeline
    For each library in --lib-list the pipeline will expect a gVCF --gvcf-dir/LIBxxxxxxx/LIBxxxxxxx.g.vcf

**Exon Mask** .bed
    To speed up the analysis of RNAseq data use the existing gene annotations and only call variants in exon regions.

Outputs
^^^^^^^^
**_SNPs[.vcf.gz,.h5]**
    Filtered single nucleotide variants
**_RefSNPs[.vcf.gz,.h5]**
    Filtered single nucleotide variants and sites called as homokaryotic reference-like
**_INDELs.vcf.gz**
    Filtered multinucleotide variants
**_SNPs_syn[.vcf.gz,.h5]**
    Synonymous sites only
**Raw.ipynb**
    Notebook reporting on QC metrics of the raw callset
**Filtered.ipynb**
    Notebook reporting on QC metrics of the filtered callset
**SNPs.ipynb**
    Notebook reporting on QC metrics of just the SNPs


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

Callset.py is designed to be run as a script taking a text file containing as list of libraries as the first arg and number of SG splits as the second.
The name of the libraries file is taken as the output-prefix for the callset.


.. code-block:: bash

    ./scripts/Callset.sh my_callset.txt 10


Tree Pipeline
-----------------

This follows on from the Callset pipeline and runs RAxML to make a phylogenetic tree.
It takes a list of libraries as an input and will create a new callset containing exactly these libraries if one does not already exist.

Key steps:
    * Apply the SNPs to the reference to create a pair of pseudohaplotype sequences
    * Extract the coding sequence regions described in the GFF
    * For each gene decide whether it has enough coverage (>80%) to be used then concatenate all genes from a sample
    * Decide if each sample has enough coverage (>80%) to be included
    * Get the 3rd codon position
    * Run RAxML and bootstrap

Flags
^^^^^
.. attribute:: --gff

    GFF file describing the gene annotations


Outputs
^^^^^^^^

.. code-block:: bash

    ./scripts/Tree.sh my_callset.txt


Transcripts Pipeline
---------------------

This uses various programs for RNAseq transcript assembly to create a set of consensus transcript annotations.

Key steps:
    * Stringtie reference guided transcript assembly
    * Cufflinks reference guided transcript assembly
    * Trinity de novo transcript assembly
    * Portcullis splice junction identification
    * Transdecoder ORF finder
    * BLAST against Uniprot
    * Mikado unification


Flags
^^^^^


Outputs
^^^^^^^^

.. code-block:: bash

    ./scripts/Transcripts.sh my_callset.txt

