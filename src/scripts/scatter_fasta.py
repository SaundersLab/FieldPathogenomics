#!/usr/bin/env python

import sys,os
import Bio.SeqIO
import numpy as np


ref_file = sys.argv[1]
output_prefix = sys.argv[2]
N_splits = int(sys.argv[3])

def scatter_fasta(ref_file, output_prefix, N_splits):
    total_length = sum([len(x) for x in Bio.SeqIO.parse(ref_file, 'fasta')])
    per_split = total_length/N_splits
    
    split_n, length_n = 0,0
    fh = open(output_prefix+'_0.interval', 'w')
    for x in Bio.SeqIO.parse(ref_file, 'fasta'):
        fh.write(x.name+"\n")
        length_n+=len(x)
        if length_n > per_split and split_n < N_splits - 1:
            fh.close()
            length_n = 0
            split_n+=1
            fh = open(output_prefix + '_'+str(split_n)+'.interval', 'w')
            
    fh.close()

if __name__ == '__main__':
    ref_file = sys.argv[1]
    output_prefix = sys.argv[2]
    N_splits = int(sys.argv[3])
    scatter_fasta(ref_file, output_prefix, N_splits)