Getting Started
=================

The code for running the fielpathogenomics pipelines is split accross three locations:

- Github https://github.com/SaundersLab/FieldPathogenomics
- Production environment (production code: /nbi/software/testing/saunderslab/FP_prod, production data: /nbi/Research-Groups/JIC/Diane-Saunders/FP_project/FP_pipeline/data, production scratch: /nbi/scratch/fieldpathogenomics)
- Local development environments

The idea being that the version of the code/data in the production environment is stable/tested and corresponds to a named (e.g. v0.2.0) version of code.
The local development environments can be used to test new changes in a way that doesn't corrupt files that other people may be using.


Setting up a local develop environment
--------------------------------------

I have written a script that automates almost all of this process.

**1. Download the init script**

    .. code-block:: bash

        ssh software
        wget https://raw.githubusercontent.com/SaundersLab/FieldPathogenomics/master/scripts/init_working_copy.sh

You will then need to specify 3 things:
    - dev_dir  The location the new develop environment will be installed
    - scratch_dir The location to store scratch files
    - prod_dir The location of the production code, this shouldn't be changed.
    - base_dir The location to put output data

**2. Run the init script**
    Use full paths please, ie no ~/

    .. code-block:: bash

        ./init_working_copy.sh --dev-dir=/usr/users/ga004/buntingd/FP_dev \
                               --base-dir=/usr/users/ga004/buntingd/FP_dev/data \
                               --scratch-dir=/nbi/scratch/buntingd \
                               --prod-dir=/nbi/software/testing/saunderslab/FP_prod 
    This will take some time

**3. Run the test scripts**

    .. code-block:: bash

        ./scripts/test.sh

    This will produce test logs which you can check and make sure everything is working.

Inside your dev_dir you will find a symlink called fieldpathogenomics that takes you to your installed version of the source code repository. You can modify this and run it in place. You can also then commit and push your changes to Github.

**4. Verify Github remote is correct**
    .. code-block:: bash

        source dev/bin/activate
        cd fieldpathogenomics
        git status
        git remote -v


Starting the Scheduler
----------------------
It is only necessary to have on instance of the scheduler running at a time, so first check if you can connect to an already running version. If not start a new one

.. code-block:: bash

    sbatch /nbi/Research-Groups/JIC/Diane-Saunders/FP_project/FP_pipeline/scripts/luigid_init.sh


This will start a new scheduler running on j256n6


Connecting to the Scheduler
---------------------------

**1. Start an SSH tunnel**

    .. code-block:: bash

        ssh -N slurm -L 8082:j256n6:8082


**2. Open your browser to http://localhost:8082**


Pulling
-------

In the scripts folder there is a script pull.sh, this will **destroy your local changes** and pull the most recent commit to Github.
A less aggressive way of merging upstream changes is through git in the normal way.

Releasing
----------
When you are done working on code in the development environment you can release it to the shared production environment using the release.sh script, this assigns a version number and creates a new, clean version of the data from your changes.

Fieldpathogenomics versions are described by a major, minor and patch number eg 0.2.4.
You should increment the patch number for small changes that do not effect existing output files, for example to fix a bug that causes the pipeline to crash.
For larger changes you need to increment the minor or major version number.
**Data is shared between versions of the code with the same patch number.**
So changing the major/minor number will cause the whole pipeline to re-run from the beginning.

.. code-block:: bash

    ./scripts/release.sh --major/--minor/--patch/--no-bump






