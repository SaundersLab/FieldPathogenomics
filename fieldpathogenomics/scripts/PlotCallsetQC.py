import pandas as pd
import sys
import re

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# def nbins(s):
#    range = s.max() - s.min()
#    n = len(s) - sum(pd.isnull(s))
#    iqr = s.quantile(0.75) - s.quantile(0.25)
#
#    return int(range * n**(1/3) /(2*iqr) )


def nbins(s):
    n = len(s) - sum(pd.isnull(s))
    return int(2 * n**(1 / 3))


def plot_info(info, df):
    df.hist(info, bins=nbins(df[info]))
    pp.savefig()


def plot_format(fmt, df):
    f = plt.figure()
    f.set_figwidth(0.5 * len(libs))
    f.set_figwidth(10)

    plt.violinplot(np.array(df[[x + '.' + fmt for x in libs]].dropna(0)))

    plt.title(fmt)
    plt.xticks(np.arange(len(libs)) + 1)
    plt.gca().set_xticklabels(libs, rotation=80)

    if fmt == 'DP':
        plt.gca().set_yscale("log")
    pp.savefig()

if __name__ == '__main__':
    info_fields = ['QD', 'DP', 'QUAL', 'FS']
    format_fields = ['GQ', 'RGQ', 'DP']

    inp, out = sys.argv[1:]
    pp = PdfPages(out)

    df = pd.read_table(inp)
    libs = [re.match('(LIB\S+)\.DP', x).groups()[0]
            for x in df.columns if re.match('(LIB\S+)\.DP', x) is not None]

    for info in info_fields:
        plot_info(info, df)

    for fmt in format_fields:
        plot_format(fmt, df)

    pp.close()
