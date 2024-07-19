from kmer_trie import KmerTrie
import pathlib
import csv
from collections import defaultdict
import os
import ast
from scipy.stats import pearsonr, linregress,kendalltau
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('csv_file', help='The name of the CSV file')
parser.add_argument('num_shards', help='total number of shards', type=int)
parser.add_argument('shard', help='shard number', type=int)
args = parser.parse_args()

nshards = int(args.num_shards)
shard = int(args.shard)

phis = [0.25, 0.5, 1, 2, 4]
methods = ['muscle', 'mafft', 'tcoffee', 'median', 'denoised', 'scaffold']

reader = csv.DictReader(open(args.csv_file))
of_dir = args.csv_file.split('.csv')[0] + '_relative/'
pathlib.Path(of_dir).mkdir(parents=True, exist_ok=True)
of = open(os.path.join(of_dir,  f'{shard}.csv'), 'w')
fields = reader.fieldnames
for phi in phis:
    fields += [m + f'_rel_{phi}_entropy' for m in methods] + [m + f'_rel_{phi}_cat_entropy' for m in methods]

writer = csv.DictWriter(of, fields)
writer.writeheader()

for i, row in enumerate(reader):
    if i % nshards != shard:
        continue
    print(i, flush=True)
    etrie = KmerTrie(12, unk_char='X')
    repeats = ast.literal_eval(row['noisy_repeats'])
    for i, r in enumerate(repeats):
        etrie.insert(r, f"r{i}")
        etrie.insert(r, "r")
    for m in methods:
        if row[m] == "":
            continue
        ent = []
        etrie.insert(row[m], m)
        for phi in phis:
            for i in range(len(repeats)):
                ent.append(etrie.relative_entropy(f"r{i}", m, phi))
            row[m + f'_rel_{phi}_entropy'] = np.mean(ent)
            row[m + f'_rel_{phi}_cat_entropy'] = etrie.relative_entropy('r', m, phi)
    writer.writerow(row)

of.close()
