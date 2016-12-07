import luigi
import os
import gzip
from collections import defaultdict
import subprocess
import pandas as pd
import hashlib
import inspect

###############################################################################
#                         Testing file emptiness                              #
###############################################################################

def isNeGz(fname):
    with gzip.open(fname, 'rb') as f:
        data = f.read(1)
    return len(data) > 0

def isNePlain(fname):
    return  os.path.getsize(fname) > 0
    
def isNeBam(fname):
    r = subprocess.run("set -o pipefail; source samtools-0.1.19 && samtools view {0} | head -c1 ".format(fname), shell=True, stdout=subprocess.PIPE)

    if r.returncode == 1:
        # samtools will return 1 if fname is not a valid bam
        # Check specifically for 1 rather than !=0 as if
        # The file is bam and non-empty head causes the retcode to be SIGPIPE 141
        raise subprocess.CalledProcessError(r.returncode, str(r.args))
        
    return len(r.stdout) > 0
    
extensions_dispath = defaultdict(lambda : isNePlain, {'bam':isNeBam,
                                                      'gz':isNeGz})
def is_not_empty(fname):

    try:
        # Try dispatching to handler determmined by file extension
        _,ext = fname.rsplit('.', 1)
        return extensions_dispath[ext](fname)

    except ValueError:
        try:
            return isNeBam(fname)
        except CalledProcessError:
            try:
                return isNeGz(fname)
            except OSError:
                return  os.path.getsize(fname) > 0
    
class CheckTargetNonEmpty(object): 
    '''This is mixin class that can be added to a luigi task to cause it to fail if the produced output is exists but is empty
        Handles checking compressed files which can have nonzero size but still be empty'''
    def complete(self):
        outputs = luigi.task.flatten(self.output())
        return super().complete() and all(map(is_not_empty, [x.path for x in outputs]))


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
        pd.Series(dict(zip(keep,[float('nan')]*len(keep))))
    
    s['Library'] = lib    
    return s

def parseLibList(liblist, base_dir):
    return pd.DataFrame([parseStarLog(os.path.join(base_dir, lib, 'Log.final.out'), lib) for lib in liblist]).set_index('Library')
    
    
###############################################################################
#                                    Git                                     #
###############################################################################
    
def file_hash(file):
    r = subprocess.run("git hash-object " + file, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    return r.stdout.strip()

def current_commit_hash(git_dir):
    r = subprocess.run("git rev-parse HEAD ", shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True, cwd=git_dir)
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
    return {f:git.git_hash(f) for f in files}
