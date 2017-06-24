function commit {

    echo "__version__ = '${NEW}'" > $version_file
    echo "Bumped version from $OLD to ${NEW}"

    ## Commit and push the version bump

ssh -t -t software << HERE
source $dev_dir/dev/bin/activate
cd $src_dir
git add $version_file
git commit -m"Version bump"
git push
git tag '${V_MAJOR}.${V_MINOR}.${V_PATCH}'
git push --tags
exit
HERE
}


## Bump the version number
cd $src_dir
#TAG=`git describe --abbrev=0 --tags`
#regex="([0-9]+)\.([0-9]+)\.([0-9]+)"

BASE_STRING=`cat $version_file`
regex="__version__ = '([0-9]+)\.([0-9]+)\.([0-9]+)'"


[[ $BASE_STRING =~ $regex ]]

V_MAJOR=${BASH_REMATCH[1]}
V_MINOR=${BASH_REMATCH[2]}
V_PATCH=${BASH_REMATCH[3]}

OLD="${V_MAJOR}.${V_MINOR}.${V_PATCH}"

if [[ $* == *--patch* ]]; then
    V_PATCH=$((V_PATCH + 1))
    NEW="${V_MAJOR}.${V_MINOR}.${V_PATCH}"
    commit

elif [[ $* == *--minor* ]]; then
    V_MINOR=$((V_MINOR + 1))
    V_PATCH=0
    NEW="${V_MAJOR}.${V_MINOR}.${V_PATCH}"
    commit

elif [[ $* == *--nobump* ]]; then
    echo "Re-releasing $OLD"
else
    echo "Please specify a release type"
    exit 1
fi


ssh -t -t software <<- HERE
    # Use the internet connected node to install required packages
    source $prod_dir/production/bin/activate
    pip install --upgrade  --force-reinstall  -e git+https://github.com/SaundersLab/FieldPathogenomics.git@${NEW}#egg=fieldpathogenomics
    exit
HERE

echo "Successfully created a new production python environment "
