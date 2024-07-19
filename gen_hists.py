import os
import nltk
import ast
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 13})
import csv
from tqdm import tqdm
import pandas as pd

import argparse
import logging

MSAS = ['denoised', 'mafft', 'muscle', 'tcoffee', 'median']
NAMES = ['SSSL', 'MAFFT', 'MUSCLE', 'T-Coffee', 'Median']

def gen_plots(fname, name, skip):
    reader = csv.DictReader(open(fname))

    edits = [defaultdict(list) for _ in MSAS]
    totals = [[] for _ in MSAS]

    for row in reader:
        no_msa = any([row[m + '_edit'] == '' for m in MSAS])
        if skip and no_msa:
            continue
        nrepeats = int(row['nrepeats'])
        for i, m in enumerate(MSAS):
            if row[m + '_edit'] == '':
                continue
            edits[i][nrepeats].append(float(row[m + '_edit']))
            totals[i].append(float(row[m + '_edit']))

    for name, total in zip(NAMES, totals):
        plt.hist(total, label=name, range=(0, 100), bins=100, histtype=u'step')
        print(name)
        print(np.average(total))
        print(np.std(total))
    plt.xlim(0, 100)
    plt.legend()
    plt.xlabel("Source Sequence Edit Distance")
    plt.ylabel("# of reads")
    plt.savefig('edit_hist.png')
    plt.clf()

    for name, edit in zip(NAMES, edits):
        print(name)
        for k in sorted(edit.keys()):
            print("{}, {}, {}".format(str(k), str(np.median(edit[k])), str(np.std(edit[k]))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='generation file', type=str, required=True)
    parser.add_argument('-n', '--name', help='data name', type=str, required=True)
    parser.add_argument('-s', '--skip', help='skip no consensus', action='store_true')
    args = parser.parse_args()
    gen_plots(args.input, args.name, args.skip)
