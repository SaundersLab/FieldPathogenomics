
# Use the internet connected node to install required packages
ssh -t -t software << HERE
source $dev_dir/dev/bin/activate
pip install --upgrade -e git+ssh://git@github.com/SaundersLab/FieldPathogenomics.git#egg=fieldpathogenomics
exit
HERE
