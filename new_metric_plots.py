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

    metrics = [defaultdict(list) for _ in METRICS]
    edits = defaultdict(list)

    for row in reader:
        if row['mafft_edit'] == '':
            continue
        base_edit = float(row['mafft_edit'])
        nrepeats = int(row['nrepeats'])
        for method in MSAS:
            if row[method + '_edit'] == '':
                continue
            #print(row[method + '_edit'])
            edits[nrepeats].append(base_edit - float(row[method + '_edit']))
            for i, metric in enumerate(METRICS):
                base_metric = float(row['mafft_' + metric])
                field = method + '_' + metric
                #print(row[field])
                metrics[i][nrepeats].append(base_metric - float(row[field]))

    for i, metric in enumerate(METRICS):
        print(metric)
        colors = cm.plasma(np.linspace(0, 1, len(metrics[i].keys())))
        for k, c in zip(sorted(metrics[i].keys()), colors):
            r = pearsonr(metrics[i][k], edits[k])[0]
            print(','.join([str(k), str(r)]))
            plt.scatter(metrics[i][k], edits[k], color=c, alpha=0.5, s=9)

        xpoints = ypoints = plt.xlim()
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('plasma'))
        sm.set_clim(vmin=2, vmax=20)
        cbar = plt.colorbar(sm)
        cbar.ax.set_yticklabels([v for v in range(2, 21, 2)])
        cbar.set_label("# Subreads")
        plt.axvline(0, color='black', linestyle='--', alpha=0.3)
        plt.axhline(0, color='black', linestyle='--', alpha=0.3)
        plt.xlabel(capwords(metric.replace('_', ' ')))
        plt.ylabel("Source Edit Distance")
        plt.savefig(metric + '.png')
        plt.clf()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='generation file', type=str, required=True)
    args = parser.parse_args()
    gen_plots(args.input)
