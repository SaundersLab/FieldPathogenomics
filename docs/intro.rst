Introduction
============

My aim is to write a set of pipelines completely covering the field pathogenomics project from reads to paper!


The Luigi Workflow Manager
**************************

Luigi is library for running complex workflows such as these. 
A luigi workflow is a graph of Task objects. The workflow is then executed by a worker under the direction of the central scheduler, we can then view the task graph and monitor execution from the web interface.
Using a luigi increases robustness and flexibility above using shell scripts etc

    * Easy to insert/remove tasks
    * Scheduler will retry tasks that fail
    * Recovery from incomplete state - each task is atomic so can resume execution after failure
    * Scheduler lazily executes the task graph - will reuse existing data if possible 
    * Easy to interface with SQL
    * It's really fast to scaffold pipelines like this, allows for rapid prototyping
    * SLURM accounting data is captured and stored allowing for a postmortem, plotting statistics about ram, cpu, disk use etc throughout the pipeline 
    
General Set-up
*******************

All the pipelines are wrapped up in a python package FieldPathogenomics and can be installed directly from Github using pip. 
This can be combined with virtualenvs to make it very simple to, starting from nothing, checkout a working version of the pipeline and run it. 
It also means it's possible to have separate development and production environments


Next Steps
***********

    In rough order of decreasing priority:
        1. Add QC/statistics to Callset
        2. Write pipelines for downstream analyses SnpEff/DAPC/RaXML/popgen/STRUCTURE
        3. Use Mikado and Portcullis to improve low qual SNPs around splice junctions <---- This would make a really nice project
        4. Comparison with other SNP calling programs and the original pipeline
        5. Phaser <- Just got published in Nature, very nice method for phasing rnaseq data
        6. Explore exciting new downstream analysis methods, I have a list
        
    Need to generally increase test and documentation coverage