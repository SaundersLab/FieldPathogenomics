import luigi
import os
import gzip

def is_not_empty(fname):
    try:
        #Raises OSError if fname has non-zero length and is not a gzip file
        with gzip.open(fname, 'rb') as f:
            data = f.read(1)
        return len(data) > 0
    except OSError:
        return  os.path.getsize(fname) > 0
    
class CheckTargetNonEmpty(object): 

    def complete(self):
        outputs = luigi.task.flatten(self.output())
        return super().complete() and all(map(is_not_empty, [x.path for x in outputs]))


