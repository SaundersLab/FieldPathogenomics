
import os
import math

import luigi
from bioluigi.slurm import SlurmExecutableTask, SlurmTask
from bioluigi.utils import CheckTargetNonEmpty

from fieldpathogenomics.utils import get_ext

import logging
logger = logging.getLogger('luigi-interface')
alloc_log = logging.getLogger('alloc_log')
alloc_log.setLevel(logging.DEBUG)

picard = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar"
python = "source /usr/users/ga004/buntingd/FP_dev/dev/bin/activate"

# Ugly hack
script_dir = os.path.join(os.path.split(__file__)[0], 'scripts')


class ScatterVCF(SlurmTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def run(self):
        # Avoid Scatter-Work-Gather-Scatter antipattern by doing a last minute check to see if
        if not self.complete():
            super().run()
        else:
            logger.info("Re-using existing scatter")

    def work(self):
        import subprocess
        import itertools
        import gzip

        zipped = get_ext(self.input().path)[1] == '.vcf.gz'
        if zipped:

            proc = subprocess.run("zgrep -ve'#' < {} |  wc -l ".format(self.input().path),
                                  shell=True, universal_newlines=True,
                                  stdout=subprocess.PIPE)

            # Create subprocesses to perform the compression
            bgzip_procs = [subprocess.Popen('/tgac/software/production/tabix/0.2.6/x86_64/bin/bgzip -c > {0}'.format(fout.path),
                                          shell=True, stdin=subprocess.PIPE, bufsize=1, universal_newlines=True)
                           for fout in self.output()]
            fouts = [p.stdin for p in bgzip_procs]

            file_context = gzip.open(self.input().path, 'rt')
        else:
            proc = subprocess.run("grep -ve'#' < {} | wc -l  ".format(self.input().path),
                                  shell=True, universal_newlines=True,
                                  stdout=subprocess.PIPE)
            fouts = [f.open('w') for f in self.output()]
            file_context = self.input().open('r')

        flen = int(proc.stdout.strip())
        print("FILE LENGTH = {}".format(flen))
        output_cycle = itertools.chain.from_iterable([itertools.repeat(p, flen // len(fouts))
                                                      for p in fouts])
        header, got_header = '', False

        with file_context as fin:
            for line in fin:
                # Catch the header
                if got_header is False:
                    if line[0] == '#':
                        header += line
                        continue
                    elif header is '':
                        raise Exception("No header on input VCF stream")
                    else:
                        got_header = True
                        for f in fouts:
                            f.write(header)
                try:
                    # Main output
                    next(output_cycle).write(line)
                except StopIteration:
                    # Clean up any remaining line to last file
                    fouts[-1].write(line)

        if zipped:
            # Create subprocesses to perform the compression
            for i, p in enumerate(bgzip_procs):
                p.stdin.flush()
                p.stdin.close()
                p.wait()
                if p.returncode != 0:
                    raise Exception("BZGIP process for file {0} failed".format(self.output()[i].path))
        else:
            for f in fouts:
                f.close()


class ScatterBED(luigi.Task, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def run(self):
        with self.input().open() as fin:
            inp = [(l, int(l.split()[2]) - int(l.split()[1])) for l in fin]

        total_seq_len = sum([x[1] for x in inp])
        perfile = math.ceil(total_seq_len / len(self.output()))

        inp_iter = iter(inp)
        files = 0
        for i, out in enumerate(self.output()):
            count = 0
            with out.open('w') as fout:
                for l, c in inp_iter:
                    fout.write(l)
                    count += c
                    if count > perfile:
                        break

                files += 1
                if files == len(self.output()):
                    # If we're on the final file dump the rest
                    fout.writelines([x[0] for x in inp_iter])


class GatherVCF(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def work_script(self):
        return '''#!/bin/bash
                picard='{picard}'
                source vcftools-0.1.13;
                source jre-8u92

                set -eo pipefail
                $picard MergeVcfs O={output}.temp.vcf.gz {in_flags}

                mv {output}.temp.vcf.gz {output}

                '''.format(picard=picard.format(mem=self.mem * self.n_cpu),
                           output=self.output().path,
                           in_flags="\\\n".join([" I= " + x.path for x in self.input()])
                           )


class GatherHD5s(SlurmTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def work(self):
        import dask.array as da
        import numpy as np
        import h5py
        from luigi.file import atomic_file

        fs = [h5py.File(f.path, mode='r') for f in self.input()]

        # Verify all H5s have the same structure
        datasets, groups, samples = [[] for x in fs], [[] for x in fs], [[] for x in fs]
        for i, f in enumerate(fs):
            f.visititems(lambda n, o: datasets[i].append(n) if isinstance(o, h5py.Dataset) else groups[i].append(n))
            samples[i] = f['samples'][:]
        if not all([set(datasets[0]) == set(x) for x in datasets]) and np.all(samples == samples[0], axis=0):
            raise Exception("All HDF5 files must have the same groups/datasets/samples!")
        datasets, groups, samples = datasets[0], groups[0], samples[0]

        # Drop Samples dataset and handle separately
        datasets = [x for x in datasets if x != 'samples']
        combined = {d: da.concatenate([da.from_array(f[d], chunks=100000) for f in fs]) for d in datasets}

        shapes = [(np.sum([f.get(d).shape for f in fs], axis=0)[0], *fs[0].get(d).shape[1:]) for d in datasets]
        dtypes = [fs[0].get(d).dtype for d in datasets]

        # Handles Samples dataset
        datasets.append('samples')
        combined.update({'samples': da.from_array(fs[0]['samples'], chunks=1)})
        shapes.append(samples.shape)
        dtypes.append(samples.dtype)

        af = atomic_file(self.output().path)
        fout = h5py.File(af.tmp_path, 'w')

        # Set up group structure
        for g in groups:
            fout.create_group(g)

        # Create the datasets
        out_datasets = {}
        for p, dtype, shape in zip(datasets, dtypes, shapes):
            g, d = os.path.split(p)
            out_datasets[p] = (fout[g] if g else fout).create_dataset(d, shape=shape, dtype=dtype,
                                                                      chunks=True, compression='gzip')
        for k in combined.keys():
            s = da.store(combined[k], out_datasets[k], compute=False)
            s.compute(num_workers=self.n_cpu)
            print("Done " + k)

        af.move_to_final_destination()


class GatherCat(SlurmExecutableTask, CheckTargetNonEmpty):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def work_script(self):
        return '''#!/bin/bash
                cat {inputs} > {output}.temp
                mv {output} {output}

                '''.format(output=self.output().path,
                           inputs=" \\\n".join([x.path for x in self.input()])
                           )


class GatherTSV(luigi.Task, CheckTargetNonEmpty):
    '''Gather TSV files making sure to only record the header line once'''

    def run(self):
        with self.output().open('w') as fout:
            for i, inp in enumerate(self.input()):
                with inp.open('r') as fin:
                    lines = fin.readlines()
                    if i == 0:
                        fout.writelines(lines)
                    else:
                        fout.writelines(lines[1:])
