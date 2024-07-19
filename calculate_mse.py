import csv
from collections import defaultdict
import os
import ast
from scipy.stats import pearsonr, linregress,kendalltau
import argparse
import pandas as pd
import numpy as np

METRICS = ['subread_edit', 'loo_edit', 'rel_0.25_entropy', 'entropy']
METHODS = ['muscle', 'tcoffee', 'mafft', 'median', 'denoised']
#METHODS = ['mafft', 'median', 'denoised']

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('csv_file', help='The name of the CSV file')
args = parser.parse_args()

data = pd.read_csv(args.csv_file, index_col='read_id')
print(data['denoised_loo_edit'])
data['denoised_loo_edit'] = data['denoised_loo_edit'].apply(lambda x: np.mean(ast.literal_eval(x)) if isinstance(x, str) else np.nan)

mses = defaultdict(dict)
pearsons = defaultdict(dict)
kendalls = defaultdict(dict)

for k in range(min(data['nrepeats']), max(data['nrepeats']) + 1):
    mask = data['nrepeats'] == k
    print("k = {}, num = {}".format(k, mask.sum()))

    for met in METRICS:
        vals = []
        edits = []
        for m in METHODS:
            vals.extend(list(data[mask][m + '_' + met]))
            edits.extend(list(data[mask][m + '_edit']))
        vals = np.asarray(vals)
        edits = np.asarray(edits)
        slope, intercept, r, p, se = linregress(vals, edits)
        preds = vals * slope + intercept
        mse = np.mean((preds - edits)**2)
        kendalls[met][k] = kendalltau(vals, edits)[0]
        mses[met][k] = mse
        pearsons[met][k] = r

for met in METRICS:
    print(met)
    #print("MSE")
    #for k in np.sort(list(mses[met].keys())):
    #    print("{}, {}".format(k, mses[met][k]))

    print("Pearson")
    for k in np.sort(list(mses[met].keys())):
        print("{}, {}".format(k, pearsons[met][k]))

    #print('\nkendall')
    #for k in np.sort(list(mses[met].keys())):
    #    print("{}, {}".format(k, kendalls[met][k]))
