sbatch --mem 24000 -N 1 -n 1 -c 1 -p nbi-long -J Callset --wrap="cd $base_dir; source $dev_dir/dev/bin/activate; \
srun python -m fieldpathogenomics.pipelines.Callset $1 $2 \
--base-dir $base_dir/data \
--workers 300"

