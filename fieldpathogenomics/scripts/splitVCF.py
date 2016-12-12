#!/usr/bin/env python

import sys,os, io
import subprocess
import itertools

base = sys.argv[1]
N_split = int(sys.argv[2])

block_size = 1000
def start_bgzip(base, i):
    proc = subprocess.Popen('/tgac/software/production/tabix/0.2.6/x86_64/bin/bgzip -c > {0}'.format(base+'_'+str(i)+".vcf.gz"),
                        shell=True, stdin=subprocess.PIPE, bufsize=1, universal_newlines=True)
    return proc

if __name__ == '__main__':
    got_header = False
    header = ''
    bgzip_procs = [start_bgzip(base, i) for i in range(N_split)]
    output_cycle = itertools.cycle(itertools.chain.from_iterable([itertools.repeat(p, block_size) for p in bgzip_procs]))

    for line in sys.stdin:
        
        # Catch the header
        if got_header is False:
            if line[0] == '#':
                header+=line
                continue
            elif header is '':
                raise Exception("No header on input VCF stream")
            else:
                got_header = True
                for p in bgzip_procs:
                    p.stdin.write(header)
                    
        next(output_cycle).stdin.write(line)
        
    for i,p in enumerate(bgzip_procs):
        p.stdin.flush()
        p.stdin.close()
        p.wait()
        if p.returncode != 0:
            raise Exception("BZGIP process for file {0} failed".format(base+'_'+str(i)+".vcf.gz"))