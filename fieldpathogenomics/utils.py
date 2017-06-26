import os
import subprocess
import pandas as pd
import hashlib
import inspect
import zlib
import logging
import time

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
gatk = "java -XX:+UseSerialGC -Xmx{mem}M -jar /tgac/software/testing/gatk/3.7.0/x86_64/GenomeAnalysisTK.jar "
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
        pd.Series(dict(zip(keep, [float('nan')] * len(keep))))

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
