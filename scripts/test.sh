
rm -rf $src_dir/tests/output;
rm -rf $src_dir/tests/scratch;

echo "Writing test logs to $dev_dir/test.{out,err}"

sbatch -n 1 --mem 4000 -p nbi-short \
-o $dev_dir/test.out \
-e $dev_dir/test.err \
--wrap="source $dev_dir/dev/bin/activate; cd $src_dir;  nosetests -vv --exe"
