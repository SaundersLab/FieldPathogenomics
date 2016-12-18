#!/usr/bin/env python

import sys
import os
import re
import io
import math
import subprocess
import pandas as pd

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

fastq_path = '/tgac/scratch/buntingd'
bam_path = '/tgac/scratch/buntingd'

pow_to_G = dict(zip(['', 'B', 'K', 'M', 'G'], [1024**-3, 1024**-3, 1024**-2, 1024**-1, 1024**0]))


def to_gigabytes(size):
    if type(size) is not str:
        return size
    base, unit = re.match("([\.\d]+)([BKMG]?)", size).groups()
    return math.ceil(float(base) * pow_to_G[unit] * 10**3) / 10**3


def parse_slurm_time(t_str):
    # format is [DD-[hh:]]mm:ss
    # return decimal hours
    s1 = t_str.split(':')
    # always have mins/secs
    secs = int(s1[-1])
    mins = int(s1[-2])

    if len(s1) > 2:
        s2 = s1[0].split('-')
        if len(s2) == 1:
            hours = int(s2[0])
            days = 0
        else:
            hours = int(s2[1])
            days = int(s2[0])
    return 24 * days + hours + mins / 60. + secs / (60.**2)


class Task():

    def __init__(self, task_id, jobid):
        super(Task).__init__()

        self.task_id = task_id
        self.type, self.lib, =  re.match('^(\S+)_(LIB\S+?)_', self.task_id).groups()
        # Differentiate Tasks run for a lib and those run for all samples
        if '_' in self.type:
            self.type = self.type.split('_')[0]
            self.lib = ''

        self.jobid = jobid

if __name__ == '__main__':
    task_file = sys.argv[1]
    base_dir = task_file.rsplit(".", 1)[0]

    tasks = {}
    with open(task_file, 'r') as f:
        for line in f:
            task_id, jobid = line.rstrip().split('\t')
            tasks[jobid] = Task(task_id, jobid)

    jobid_query = ','.join([t.jobid + '.0' for t in tasks.values()])
    p = subprocess.run("sacct -P --format=jobid,elapsed,MaxDiskWrite,MaxDiskRead,AveRSS,MaxRSS,AveVMSize,MaxVMSize,state,ExitCode -j " +
                       jobid_query, shell=True, universal_newlines=True, stdout=subprocess.PIPE)
    task_table = pd.read_table(io.StringIO(p.stdout), sep='|')

    task_table['JobID'] = task_table['JobID'].astype(str).str.split('.').str.get(0)
    task_table['Lib'] = [tasks[jid].lib for jid in task_table['JobID']]
    task_table['Type'] = [tasks[jid].type for jid in task_table['JobID']]

    task_table['Elapsed'] = task_table['Elapsed'].map(parse_slurm_time)
    for mem in ['MaxDiskWrite', 'MaxDiskRead', 'AveRSS', 'MaxRSS', 'AveVMSize', 'MaxVMSize']:
        task_table[mem] = task_table[mem].map(to_gigabytes)

    completed = task_table.ix[task_table['State'] == 'COMPLETED']

    # By Res
    res_path = os.path.join(base_dir, 'by_res')
    os.makedirs(res_path, exist_ok=True)

    for res in ['MaxDiskWrite', 'MaxDiskRead', 'AveRSS', 'MaxRSS', 'AveVMSize', 'MaxVMSize']:
        completed.boxplot(column=res, by='Type', )
        plt.xticks(rotation=75)
        plt.ylabel('Gigabytes')
        plt.savefig(os.path.join(res_path, res + ".pdf"))
        plt.close()

    completed.boxplot(column='Elapsed', by='Type', )
    plt.xticks(rotation=75)
    plt.ylabel('Hours')
    plt.savefig(os.path.join(res_path, "Elapsed.pdf"))
    plt.close()

    # By Task
    task_path = os.path.join(base_dir, 'by_task')
    os.makedirs(task_path, exist_ok=True)

    for task in set(completed['Type']):
        completed[completed['Type'] == task].boxplot(
            column=['MaxDiskWrite', 'MaxDiskRead', 'AveRSS', 'MaxRSS', 'AveVMSize', 'MaxVMSize'])
        plt.xticks(rotation=75)
        plt.title(task)
        plt.ylabel('Gigabytes')
        plt.savefig(os.path.join(task_path, task + ".pdf"))
        plt.close()

    # Elapsed time breakdown
    completed.groupby(['Type', 'Lib'])['Elapsed'].sum().unstack(
        'Type').plot(kind='bar', stacked=True, figsize=(20, 10))
    plt.savefig(os.path.join(base_dir, 'Lib_time.pdf'))

    # By fastq size
    fastq_size = []
    os.makedirs(os.path.join(base_dir, 'scatter_plots'), exist_ok=True)
    for lib in completed['Lib']:
        if lib:
            fastq_size.append(os.path.getsize(os.path.join(
                fastq_path, lib, 'raw_R1.fastq.gz')) / 1024**3)
        else:
            fastq_size.append(None)

    completed['Fastq_size'] = fastq_size

    gb = completed.groupby('Type')

    pp = PdfPages(os.path.join(base_dir, 'scatter_plots', 'fastq.pdf'))
    plt.figure()
    for name, group in gb:
        if name in ['FetchFastqGZ', 'PythonFilter', 'FastxTrimmer', 'Star']:
            plt.plot(group.Fastq_size, group.Elapsed, '.', label=name)
    plt.xlabel('Fastq.gz size (GB)')
    plt.ylabel('Time (Hours)')
    plt.legend()
    pp.savefig()

    plt.figure()
    for name, group in gb:
        if name in ['FetchFastqGZ', 'PythonFilter', 'FastxTrimmer', 'Star']:
            plt.plot(group.Fastq_size, group.AveRSS, '.', label=name)
    plt.xlabel('Fastq.gz size (GB)')
    plt.ylabel('AveRSS (GB)')
    plt.legend()
    pp.savefig()

    plt.figure()
    for name, group in gb:
        if name in ['FetchFastqGZ', 'PythonFilter', 'FastxTrimmer', 'Star']:
            plt.plot(group.Fastq_size, group.AveVMSize, '.', label=name)
    plt.xlabel('Fastq.gz size (GB)')
    plt.ylabel('AveVMSize (GB)')
    plt.legend()
    pp.savefig()
    pp.close()

    # By bam size
    bam_size = []
    for lib in completed['Lib']:
        if lib:
            bam_size.append(os.path.getsize(os.path.join(
                bam_path, lib, 'Aligned.out_cleaned.bam')) / 1024**3)
        else:
            bam_size.append(None)

    completed['Bam_size'] = bam_size

    pp = PdfPages(os.path.join(base_dir, 'scatter_plots', 'bam.pdf'))

    plt.figure()
    for name, group in gb:
        if name in ['AddReadGroups', 'CleanSam', 'HaplotypeCaller', 'MarkDuplicates', 'PlotAlleleFreq', 'SplitNCigarReads']:
            plt.plot(group.Fastq_size, group.Elapsed, 'o', label=name)
    plt.xlabel('BAM size (GB)')
    plt.ylabel('Time (Hours)')
    plt.legend()
    pp.savefig()

    plt.figure()
    for name, group in gb:
        if name in ['AddReadGroups', 'CleanSam', 'HaplotypeCaller', 'MarkDuplicates', 'PlotAlleleFreq', 'SplitNCigarReads']:
            plt.plot(group.Fastq_size, group.AveRSS, 'o', label=name)
    plt.xlabel('BAM size (GB)')
    plt.ylabel('AveRSS (GB)')
    plt.legend()
    pp.savefig()

    plt.figure()
    for name, group in gb:
        if name in ['AddReadGroups', 'CleanSam', 'HaplotypeCaller', 'MarkDuplicates', 'PlotAlleleFreq', 'SplitNCigarReads']:
            plt.plot(group.Fastq_size, group.AveVMSize, 'o', label=name)
    plt.xlabel('BAM size (GB)')
    plt.ylabel('AveVMSize (GB)')
    plt.legend()
    pp.savefig()

    pp.close()
