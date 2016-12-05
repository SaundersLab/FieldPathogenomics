Protocols
===========

Adding new samples
-------------------

1. Add the metadata to the database
+++++++++++++++++++++++++++++++++++
TODO

2. Generate a text file containing the library names
+++++++++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: bash

    echo 'LIB22849
          LIB23170
             .
             .
             .
          LIB23262' > libs.txt
          
3. Run the Library pipeline
++++++++++++++++++++++++++++
Example slurm wrapper script

.. code-block:: bash

   #!/bin/sh
   #SBATCH --mem 20000 
   #SBATCH -n 1
   #SBATCH -N 1
    
    cd /tgac/workarea/collaborators/saunderslab/FP_pipeline
    source production/bin/activate
    export LUIGI_CONFIG_PATH=/tgac/workarea/collaborators/saunderslab/FP_pipeline/luigi.cfg

    srun python -m src.pipelines.Library libs.txt \
    --base-dir /tgac/workarea/collaborators/saunderslab/FP_pipeline/data/ \
    --scratch-dir /tgac/scratch/buntingd/ \
    --read-dir /tgac/data/read/*DianeSaunders* \
    --workers 250 

4. Verify QC results and alignments stats
++++++++++++++++++++++++++++++++++++++++++


Generating a new callset
-------------------------

1. Run the Library pipeline 
+++++++++++++++++++++++++++++
Need to have a gVCF for every library and have them all in the same folder

2. Make a list of the libraries to include
+++++++++++++++++++++++++++++++++++++++++++

3. Make a BED file of the regions to genotype
+++++++++++++++++++++++++++++++++++++++++++++

4. Run the Callset pipeline
+++++++++++++++++++++++++++++

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
    
    
5. Visualise the callset in IGV
+++++++++++++++++++++++++++++++
Make sure it looks sensible

6. Check QC plots and stats
++++++++++++++++++++++++++++
TODO


Creating a new fieldpathogenomics environment
----------------------------------------------

Please use python virtual environments to keep the version of packages consistent!!!

1. Create a new python virtual environment
++++++++++++++++++++++++++++++++++++++++++
Starting in the folder you want to create the env in:
Call the new virtualenv prod

.. code-block:: bash
    $ source python-3.5.1
    $ python -m virtualenv -p `which python3` prod
2. Modfiy the virtualenv wrapper
++++++++++++++++++++++++++++++++++++
This is just for convenience as on the cluster we have to use source.
Add the line `source python-3.5.1; source git-1.8.1.2; export TMPDIR=/tgac/scratch/buntingd/` to the top of the activate script `prod/bin/activate`

3. Install FieldPathogenomics
++++++++++++++++++++++++++++++

pip install git+https://github.com/SaundersLab/FieldPathogenomics.git

Optionally install a specific commit/branch

