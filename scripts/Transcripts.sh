sbatch --mem 8000 -N 1 -n 1 -c 1 -p nbi-long -J Transcripts --wrap="cd $base_dir; source $dev_dir/dev/bin/activate; \
srun python -m fieldpathogenomics.pipelines.Transcripts $1 \
--base-dir $base_dir \
--scratch-dir $scratch_dir \
--workers 300"

