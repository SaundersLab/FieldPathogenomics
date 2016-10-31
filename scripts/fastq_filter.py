#!/usr/bin/env python

import Bio.SeqIO
import gzip
import argparse


class no_Ns():
    '''Filter any read that contains an N'''
    def __init__(self):
        self.failures = 0
        
    def __call__(self, r1, r2):
        predicate = not (('N' in r1.seq) or ('N' in r2.seq))
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
        self.L = L
        self.failures = 0
    def __call__(self, r1, r2):
        predicate = (len(r1.seq) == self.L) and (len(r2.seq) == self.L)         
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
        with gzip.open(self.in_r1, 'rt') as r1_in_fh, gzip.open(self.in_r2, 'rt') as r2_in_fh:
            with gzip.open(self.out_r1, 'wt') as r1_out_fh, gzip.open(self.out_r2, 'wt') as r2_out_fh:
                in_count, out_count = 0,0
                for r1, r2 in zip(Bio.SeqIO.parse(r1_in_fh, 'fastq'), Bio.SeqIO.parse(r2_in_fh, 'fastq')): 
                    in_count += 1
                    if all([f(r1, r2) for f in self.filters]):
                        out_count += 1
                        Bio.SeqIO.write(r1, r1_out_fh, 'fastq')
                        Bio.SeqIO.write(r2, r2_out_fh, 'fastq')
                         
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
    
