import torch
import sys
import os
from tqdm import tqdm
import multiprocessing as mp
import nltk
import csv
import ast
from copy import copy
import json
import matplotlib.pyplot as plt
import numpy as np
from denoising.data_utils import decode, codon_tok
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict
from scipy.stats import ttest_1samp
from kmer_trie import KmerTrie
from itertools import product
import argparse

def kmer_entropy(seqs, k, phis, compare_seqs, all_k):
    res = []
    kmer_counts = KmerTrie(k)

    # build kmer dict
    for seq in seqs:
        for subseq in seq.split():
            kmer_counts.insert(subseq)

    # get entropy
    for cseq in compare_seqs:
        if cseq == '':
            res.append(None)
        else:
            kmer_counts.insert(cseq)
            res.append(kmer_counts.get_entropy(phis, all_k))
            kmer_counts.delete(cseq)

    return res

def row_worker(fin):
    rows, start_id, k, phis, all_k, header = fin
    res = []
    print("Starting process with {} rows...".format(len(rows)))
    ridx = header.index('read')
    sidx = header.index('scaffold')
    didx = header.index('denoised')
    cidx = header.index('consensus')
    for i, row in enumerate(rows):
        print(start_id + i, flush=True)
        repeats = ast.literal_eval(row[ridx])
        scaffold = row[sidx]
        denoised = row[didx]
        consensus = row[cidx]
        tsent, tdent, tcent = kmer_entropy(repeats, k, phis, [scaffold, denoised, consensus], all_k)
        res.append([start_id + i, len(repeats), tsent, tcent, tdent, row[-2], row[3]])
    return res

def gen_entropy(fname, k, phis, pool_size, all_k):
    pool = mp.Pool(pool_size)

    reader = csv.reader(open(fname, 'r'))
    header = next(reader)
    rows = list(reader)
    nrows = len(rows)

    chunks = nrows // pool_size + 1
    fins = [(rows[i * chunks:(i + 1) * chunks], i*chunks, k, phis, all_k, header) for i in range(pool_size)]

    res = []
    for r in pool.map(row_worker, fins):
        res.extend(r)

    outf = os.path.join(os.path.dirname(fname), 'res_kmer_{}_{}.json'.format(str(k), '_'.join([str(p) for p in phis])))
    json.dump(res, open(outf, 'w'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='input generation file')
    parser.add_argument('-k', help='max kmer limit', type=int)
    parser.add_argument('-p', '--phi', nargs='+', help='phi parameter', type=float)
    parser.add_argument('-s', '--size', help='pool size', type=int)
    parser.add_argument('-a', '--all', help='all ks', action='store_true')
    args = parser.parse_args()
    gen_entropy(args.file, args.k, args.phi, args.size, args.all)
