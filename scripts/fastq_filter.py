#!/usr/bin/env python

import gzip
import argparse
import io


class no_Ns():
    '''Filter any read that contains an N'''
    def __init__(self):
        self.failures = 0
        
    def __call__(self, r1, r2):
        predicate = not (('N' in r1) or ('N' in r2))
        if not predicate:
            self.failures += 1
        return predicate

    def __str__(self):
        return "No N's: {0}".format(self.failures)
        
class exact_length():
    '''Filter any read that does not have length ``L``
    :param L: Read filter length
    '''
    def __init__(self, L=101):
        self.L = L + 1 #Add 1 to count the newline
        self.failures = 0
        
    def __call__(self, r1, r2):
        predicate = (len(r1) == self.L) and (len(r2) == self.L)         
        if not predicate:
            self.failures += 1
        return predicate
        
    def __str__(self):
        return "Exact length: {0}".format(self.failures)
        
class FastqFilter():
    def __init__(self, in_r1, in_r2, out_r1, out_r2, filters):
        
        self.in_r1 = in_r1
        self.in_r2 = in_r2
        self.out_r1 = out_r1
        self.out_r2 = out_r2
        self.filters = filters
        
        self.apply()

    def apply(self):
        
        ## Use mid sized caches 134mb to reduce network load
        with open(self.in_r1, 'rb', buffering=2**27) as r1_in_gz, open(self.in_r2, 'rb', buffering=2**27) as r2_in_gz:
            r1_in_fh = io.TextIOWrapper(gzip.GzipFile(mode='r', fileobj=r1_in_gz))
            r2_in_fh = io.TextIOWrapper(gzip.GzipFile(mode='r', fileobj=r2_in_gz))
            
            ## Fully buffer the output files, so script is (more) atomic
            r1_buf, r2_buf = io.BytesIO(), io.BytesIO()
            r1_out_gz, r2_out_gz = gzip.GzipFile(mode='w', fileobj=r1_buf), gzip.GzipFile(mode='w', fileobj=r2_buf)
            r1_out_text, r2_out_text = io.TextIOWrapper(r1_out_gz),io.TextIOWrapper(r2_out_gz)
            
            in_count, out_count = 0,0
            zipped = zip(r1_in_fh, r2_in_fh)
            for r1_header, r2_header in zipped:
                 
                if r1_header[0] != '@' or r2_header[0] != '@':
                    print("Not at the start of a record!!")
                    continue
                
                r1_seq, r2_seq = next(zipped)
                r1_desc, r2_desc = next(zipped)
                r1_qual, r2_qual = next(zipped)
                
                in_count += 1
                if all([f(r1_seq, r2_seq) for f in self.filters]):
                    out_count += 1
                    r1_out_text.writelines([r1_header, r1_seq, r1_desc, r1_qual])
                    r2_out_text.writelines([r1_header, r1_seq, r1_desc, r1_qual])
                
                if in_count > 1000:
                    break
                    
                if in_count % 100000 == 0 :
                    print("Done {0}".format(in_count))
                    
            with open(self.out_r1, 'wb') as r1_out_fh, open(self.out_r2, 'wb',) as r2_out_fh:
                r1_out_gz.close()
                r2_out_gz.close()
                r1_out_fh.write(r1_buf.getvalue())
                r2_out_fh.write(r2_buf.getvalue())
                
        print("In : {0} Failed: ".format(in_count) + "\t".join([str(f) for f in self.filters]))

    
if __name__ == '__main__':
        
    parser = argparse.ArgumentParser()
    parser.add_argument('in_R1')
    parser.add_argument('in_R2')
    parser.add_argument('out_R1')
    parser.add_argument('out_R2')
    parser.add_argument('-L', required=False, default=101)
    args = parser.parse_args()
    
    FastqFilter(args.in_R1, args.in_R2, args.out_R1, args.out_R2, [no_Ns(), exact_length(int(args.L))])
