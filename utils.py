import luigi
import os
import gzip
from collections import defaultdict
import subprocess

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
        raise subprocess.CalledProcessError()
        
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

    def complete(self):
        outputs = luigi.task.flatten(self.output())
        return super().complete() and all(map(is_not_empty, [x.path for x in outputs]))


