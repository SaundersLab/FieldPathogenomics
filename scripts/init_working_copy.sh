#!/bin/bash -e

scratch_dir=/tgac/scratch/buntingd/   # <------------- Modify this
dev_dir=/usr/users/ga004/buntingd/test # <------------- Modify this
prod_dir=/usr/users/ga004/buntingd/FP_prod   #/nbi/Research-Groups/JIC/Diane-Saunders/FP_pipeline
version_file=$prod_dir/production/src/fieldpathogenomics/fieldpathogenomics/version.py


function venv_create {

    activate_prefix='source python-3.5.1;
    source git-1.8.1.2;
    export TMPDIR=$scratch_dir'

    # Create the new virtualenv
    mkdir -p $dev_dir
    cd $dev_dir
    source python-3.5.1;
    virtualenv -p `which python3` dev

    # Add the source commands and environment variable defs to the activate script
    echo $activate_prefix > temp
    cat dev/bin/activate >> temp
    mv -f temp dev/bin/activate
}


function install_fieldpathogenomics {
    # Pull the code from the tip of master and install

ssh -t -t software << HERE
# Use the internet connected node to install required packages
source $dev_dir/dev/bin/activate
pip install --upgrade  --force-reinstall  -e git+https://github.com/SaundersLab/FieldPathogenomics.git@master#egg=fieldpathogenomics
exit
HERE

    echo "Installed fieldpathogenomics to $dev_dir/dev/src/fieldpathogenomics/fieldpathogenomics"

}

function install_requirements {
    # Use the requirements.txt to install python packages
ssh -t -t software << HERE
# Use the internet connected node to install required packages
source $dev_dir/dev/bin/activate
pip install -r $dev_dir/dev/src/fieldpathogenomics/requirements.txt
exit
HERE

    echo "Installed the requirements in  $dev_dir/dev/src/fieldpathogenomics/requirements.txt"

}

function install_scripts {
    # Copy and localise the supporting scripts
    # This makes the following vars available to all scripts:
    # $dev_dir
    # $src_dir
    # $prod_dir
    # $scratch_dir

    cp -fr $dev_dir/dev/src/fieldpathogenomics/scripts $dev_dir/scripts
    cd $dev_dir/scripts

vars=$(printf "#!/bin/bash -e
prod_dir=$prod_dir
dev_dir=$dev_dir
scratch_dir=$scratch_dir
src_dir=$dev_dir/dev/src/fieldpathogenomics/")

    printf '%b\n' "$vars" | cat - release.sh > temp && mv temp release.sh
    printf '%b\n' "$vars" | cat - test.sh > temp && mv temp test.sh
    printf '%b\n' "$vars" | cat - pull.sh > temp && mv temp pull.sh

}

venv_create;

install_fieldpathogenomics;

install_requirements;

install_scripts
