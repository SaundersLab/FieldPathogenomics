sbatch --mem 8000 -N 1 -n 1 -c 1 -p nbi-long -J Callset --wrap="cd $base_dir; source $dev_dir/dev/bin/activate; \
srun python -m fieldpathogenomics.pipelines.Tree $1 \
--consensus-type $2 \
--base-dir $base_dir/data \
--scratch-dir $scratch_dir \
--workers 300"
