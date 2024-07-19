import csv
from string import capwords
from matplotlib import ticker
import matplotlib.cm as cm
from collections import defaultdict
import ast
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
from scipy.stats import pearsonr, ttest_rel

MSAS = ['muscle', 'tcoffee', 'median', 'quickmedian']
METRICS = ['loo_edit', 'entropy', 'subread_edit']

def gen_plots(fname):
    reader = csv.DictReader(open(fname))

    corrs = [defaultdict(list) for _ in METRICS]
    edits = defaultdict(list)

    for row in reader:
        if row['mafft_edit'] == '':
            continue
        nrepeats = int(row['nrepeats'])
        for i, metric in enumerate(METRICS):
            edits = []
            metrics = []
            for method in MSAS:
                if row[method + '_edit'] == '' or row[method + '_' + metric] == '':
                    continue
                edits.append(float(row[method + '_edit']))
                metrics.append(float(row[method + '_' + metric]))
            corrs[i][nrepeats].append(pearsonr(edits, metrics)[0]**2)

    for i, metric in enumerate(METRICS):
        print(metric)
        for k in sorted(corrs[i].keys()):
            print(','.join([str(k), str(np.mean(corrs[i][k]))]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='generation file', type=str, required=True)
    args = parser.parse_args()
    gen_plots(args.input)
