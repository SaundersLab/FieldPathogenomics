#!/bin/bash -e

scratch_dir=/tgac/scratch/buntingd/
install_dir=/usr/users/ga004/buntingd/test
env_name="dev"


function venv_create {

    activate_prefix='source python-3.5.1;
    source git-1.8.1.2;
    export TMPDIR=$scratch_dir'

    # Create the new virtualenv
    mkdir -p $install_dir
    cd $install_dir
    source python-3.5.1;
    virtualenv -p `which python3` $env_name

    # Add the source commands and environment variable defs to the activate script
    echo $activate_prefix > temp
    cat $env_name/bin/activate >> temp
    mv -f temp $env_name/bin/activate
}


function install_fieldpathogenomics {
    # Pull the code from the tip of master and install

ssh -t -t software << HERE
# Use the internet connected node to install required packages
source $install_dir/$env_name/bin/activate
pip install --upgrade  --force-reinstall  -e git+https://github.com/SaundersLab/FieldPathogenomics.git@master#egg=fieldpathogenomics
exit
HERE

    echo "Installed fieldpathogenomics to $install_dir/$env_name/src/fieldpathogenomics/fieldpathogenomics"

}

function install_requirements {
    # Use the requirements.txt to install python packages
ssh -t -t software << HERE
# Use the internet connected node to install required packages
source $install_dir/$env_name/bin/activate
pip install -r $install_dir/$env_name/src/fieldpathogenomics/requirements.txt
exit
HERE

    echo "Installed the requirements in  $install_dir/$env_name/src/fieldpathogenomics/requirements.txt"

}

function install_scripts {
    # Copy and localise the supporting scripts
    cp -fr $install_dir/$env_name/src/fieldpathogenomics/scripts $install_dir/scripts
    cd $install_dir/scripts

    vars=$(printf "#!/bin/bash -e\ndev_dir=$install_dir\nsrc_dir=$install_dir/$env_name/src/fieldpathogenomics/\n")
    echo "$vars" | cat - release.sh > temp && mv temp release.sh

}

venv_create;

install_fieldpathogenomics;

install_requirements;

install_scripts
