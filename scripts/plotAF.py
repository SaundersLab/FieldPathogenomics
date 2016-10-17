#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys, re
import matplotlib 
matplotlib.use('pdf')
import matplotlib.pyplot as plt

input_tsv = sys.argv[1]
output = sys.argv[2]

# Tab separated file with columns 'CHROM', 'POS', 'REF', 'ALT', 'DP', 'LIBxxxx.AD'

tsv = pd.read_table(input_tsv, sep='\t')

# If multiple libraries are present, just pull out the first
ad_col_name = tsv.filter(regex='\S+AD').columns[0]
lib_name = re.match("(\S+)\.AD", ad_col_name).groups()[0]

counts = [np.sort([float(y) for y in x.split(',')])[::-1] for x in tsv[ad_col_name]]
freqs = [x/np.sum(x) for x in counts]
df = pd.DataFrame(freqs)
df[df==0] = float('nan')

ax_list = df.hist(sharex=True, sharey=True, range=(0,1), bins=20)
fig = ax_list.flat[0].get_figure()
plt.gcf().suptitle(lib_name,  fontsize=20)
plt.gcf().text(0.5, 0.04, 'Allele frequency', ha='center')
plt.gcf().text(0.02, 0.5, 'Counts', va='center', rotation='vertical')

plt.gcf().savefig(output)
