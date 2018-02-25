import os
import subprocess
import pandas as pd
import hashlib
import inspect
import zlib
import logging
import time
import allel
import numpy as np
import Bio
import Bio.SeqIO
import Bio.AlignIO
import itertools
from collections import Counter

###############################################################################
#                           Variant Indexing                                  #
###############################################################################


class Contig():
    '''Makes Contigs follow the correct sort order by skipping a prefix
       of length plen'''

    def __init__(self, contig, plen):
        self.plen = plen
        self.contig = contig

    def __delitem__(self, key):
        self.contig.__delitem__(key)

    def __getitem__(self, key):
        return self.contig.__getitem__(key)

    def __setitem__(self, key, value):
        self.contig.__setitem__(key, value)

    def __eq__(self, *args, **kwargs):
        return super().__eq__(*args, **kwargs)

    def __lt__(self, other):
        return int(self[self.plen:]).__lt__(int(other[self.plen:]))

    def __le__(self, other):
        return int(self[self.plen:]).__le__(int(other[self.plen:]))

    def __gt__(self, other):
        return int(self[self.plen:]).__gt__(int(other[self.plen:]))

    def __ge__(self, other):
        return int(self[self.plen:]).__ge__(int(other[self.plen:]))

    def __repr__(self):
        return self.contig.__repr__()


def index_variants(variants, plen):
    '''Create a CHROM/POS multiindex for VariantsTable
       where the contig number has a prefix (eg PST130_)
       of length plen'''
    index = allel.SortedMultiIndex(np.array([Contig(x, plen) for x in variants['CHROM'][:]], dtype=Contig), variants['POS'][:])
    return index


iupac_het = {'AG': 'R', 'GA': 'R', 'AT': 'W', 'TA': 'W',
             'CT': 'Y', 'TC': 'Y', 'GT': 'K', 'TG': 'K',
             'GC': 'S', 'CG': 'S', 'AC': 'M', 'CA': 'M'}


def memoized(f):
    cache = {}
    def ret(*args):
        if args not in cache:
            cache[args] = f(*args)
        else:
            cache[args], r = itertools.tee(cache[args])
        return r

    return ret


class GeneSet():
    def __init__(self, genes, genotypes, variants, samples, index):

        self._genotypes = genotypes
        self._variants = variants

        self.samples = samples
        self.index = index

        self.masks, gene_idxs, self.lengths = [], [], []
        for i, c in enumerate(((g[0], g[3], g[4]) for g in genes)):
            try:
                self.masks.append(index.locate_range(*c))
                gene_idxs.append(i)
                self.lengths.append(c[2] - c[1] + 1)
            except KeyError:
                pass
        # Make sure we only retain genes that have some
        # sites associated with them to avoid KeyErrors
        # downstream
        self.genes = genes[gene_idxs]

    def filter(self, filter):
        return GeneSet(self.genes[list(filter)], self._genotypes, self._variants, self.samples, self.index)

    def nucleotide_diversity(self, pop_samples=None):
        subpop = [self.samples.index(x) for x in pop_samples] if pop_samples else None

        for i, (gt, l) in enumerate(zip(self.genotypes(), self.lengths)):
            mpd = allel.mean_pairwise_difference(gt.count_alleles(subpop=subpop)).sum()
            yield mpd / l

    def genotypes(self):
        return (self._genotypes[m] for m in self.masks)

    def variants(self):
        return (self._variants[m] for m in self.masks)

    def site_coverages(self):
        return (g.is_called().sum(axis=1) / len(self.samples)for g in self.genotypes())

    @memoized
    def allele_counts(self, subpop=None):
        if subpop and isinstance(subpop[0], str):
            subpop = [self.samples.index(x) for x in subpop]
        return (gt.count_alleles(subpop=subpop) for gt in self.genotypes())

    @memoized
    def sample_coverages(self):
        return (g.is_called().sum(axis=0) / g.shape[0] for g in self.genotypes())

    def invariants(self):
        for var in self.variants():
            yield {k: v for k, v in
                   Counter(var['REF'][~var['is_snp']]).items()
                   if len(k) == 1}

    def variant_genotypes(self):
        for g, v in zip(self.genotypes(), self.variants()):
            # Dask boolean indexing seems to break if the dimension is
            # 1 because it compares True >= dim so hack around that
            if g.shape[0] == 1 and v['is_snp'][0]:
                yield g
            else:
                yield g[v['is_snp']]

    def MSA(self):
        for v, g, vg in zip(self.variants(), self.genotypes(), self.variant_genotypes()):

            msa = [Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(''), id=s)
                   for s in self.samples]
            alts = v[v['is_snp']]['ALT'][:, 0]
            refs = v[v['is_snp']]['REF']
            Ns = np.array(['N'] * vg.shape[0])
            hets = np.array([iupac_het[x[0] + x[1]] for x in zip(refs, alts)])
            X = np.vstack((alts, hets, refs, Ns)).T

            for i in range(vg.shape[0]):
                site = X[i, vg.to_n_ref(fill=-1)[i]]

                # RAxML counts a site with just N and a base as invariant so remove these
                uniq_bases = set(list(site))
                if len(uniq_bases) == 1 or ('N' in uniq_bases and len(uniq_bases) == 2):
                    continue

                for j in range(vg.shape[1]):
                    msa[j] += site[j]
            yield Bio.Align.MultipleSeqAlignment(msa)



###############################################################################
#                               File handling                                  #
###############################################################################


def get_ext(path):
    '''Split path into base and extention, gracefully handling compressed extensions eg .gz'''
    base, ext1 = os.path.splitext(path)
    if ext1 == '.gz':
        base, ext2 = os.path.splitext(base)
        return base, ext2 + ext1
    else:
        return base, ext1

###############################################################################
#                               Python paths                                   #
###############################################################################


python = "source " + os.environ['VIRTUAL_ENV'] + "/bin/activate"
notebooks = os.path.join(os.path.split(__file__)[0], 'notebooks')
reference_dir = '/nbi/Research-Groups/JIC/Diane-Saunders/FP_project/FP_pipeline/reference'

###############################################################################
#                               Java paths                                    #
###############################################################################

picard = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar"
gatk = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/gatk/3.8.0/x86_64/GenomeAnalysisTK.jar "
trimmomatic = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/trimmomatic/0.36/x86_64/bin/trimmomatic-0.36.jar "
snpeff = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/snpeff/4.3g/x86_64/snpEff.jar "
snpsift = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/snpeff/4.3g/x86_64/SnpSift.jar "

###############################################################################
#                             Pipeline init                                   #
###############################################################################


def logging_init(log_dir, pipeline_name):
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger('luigi-interface')
    alloc_log = logging.getLogger('alloc_log')

    logging.disable(logging.DEBUG)
    timestr = time.strftime("%Y%m%d-%H%M%S")

    fh = logging.FileHandler(os.path.join(log_dir, pipeline_name + "_" + timestr + ".log"))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    alloc_fh = logging.FileHandler(os.path.join(log_dir, pipeline_name + "_" + timestr + ".salloc.log"))
    alloc_fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    alloc_fh.setFormatter(formatter)
    alloc_log.addHandler(alloc_fh)

    return logger, alloc_log


###############################################################################
#                                 Checksum                                   #
###############################################################################


def checksum(file):
    '''Adler32 checksum. Uses blocksize of 500MB. Takes ~30s to checksum 10GB'''
    BLOCK_SIZE = 2**29
    with open(file, 'rb') as f:
        value = 0
        while True:
            data = f.read(BLOCK_SIZE)
            if not data:
                break
            value = zlib.adler32(data, value)
    return value


###############################################################################
#                           Parsing STAR logs                                  #
###############################################################################
keep = [
    'Number of input reads',
    'Average input read length',
    'Uniquely mapped reads number',
    'Uniquely mapped reads %',
    'Average mapped length',
    'Mismatch rate per base, %',
    'Started job on',
]


def parseStarLog(logfile, lib):
    '''Parse the star Log.final.out file :param: logfile for the library :param: lib
       returns a pandas.Series of the fields defined in keep'''
    try:
        s = pd.Series()
        with open(logfile, 'r') as f:
            for line in f:
                split = line.strip().split(" |\t")
                if len(split) == 2 and split[0] in keep:
                    s[split[0]] = float(split[1][:-1]) if '%' in split[1] else split[1]
    except FileNotFoundError:
        s = pd.Series(dict(zip(keep, [float('nan')] * len(keep))))

    s['Library'] = lib
    return s


def parseLibList(liblist, base_dir):
    return pd.DataFrame([parseStarLog(os.path.join(base_dir, lib, 'Log.final.out'), lib) for lib in liblist]).set_index('Library')


###############################################################################
#                                    Git                                     #
###############################################################################

def file_hash(file):
    r = subprocess.run("git hash-object " + file, shell=True, check=True,
                       stdout=subprocess.PIPE, universal_newlines=True)
    return r.stdout.strip()


def current_commit_hash(git_dir):
    r = subprocess.run("git rev-parse HEAD ", shell=True, check=True,
                       stdout=subprocess.PIPE, universal_newlines=True, cwd=git_dir)
    return r.stdout.strip()

###############################################################################
#                          Pipeline hashing                                   #
###############################################################################


def ancestor(task):
    '''Return a set of all Tasks that are ancestors of :param: task'''
    tree = set()
    if task.deps() != set():
        for d in task.deps():
            tree.add(d)
            tree.update(ancestor(d))
    return list(tree)


def hash_pipeline(task):
    '''Creates a hash of the source code of :param: task and all tasks upstream of it'''
    sha = hashlib.sha1()
    task_list = [task] + ancestor(task)
    sources = sorted([inspect.getsource(type(t)) for t in task_list])
    for s in sources:
        sha.update(s.encode())

    return sha.hexdigest()


def hash_files(task):
    '''Get all of the files used in any task upstream of :parm: task and returns there git file hash'''
    task_list = [task] + list(ancestor(task))
    files = [inspect.getfile(type(t)) for t in task_list]
    return {f: git.git_hash(f) for f in files}
