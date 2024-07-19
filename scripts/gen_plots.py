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

import argparse
import logging
logger = logging.getLogger('align')

def gen_plots(fname, name):
    reader = csv.reader(open(fname))
    header = next(reader)
    consensus_idx = header.index('consensus')
    scaffold_idx = header.index('scaffold')

    dedit_idx = header.index('denoised_r_edit')
    cedit_idx = header.index('consensus_r_edit')
    dmatch_idx = header.index('denoised_match')
    cmatch_idx = header.index('consensus_match')
    r_idx = header.index('read')

    nrepeats_idx = defaultdict(lambda: defaultdict(list))
    totals = defaultdict(list)

    ns = []
    lens = []

    for row in reader:
        repeats = eval(row[r_idx])
        n_repeats = len(repeats)
        ns.append(n_repeats)
        lens.extend([len(r) for r in repeats])
        n_repeats = len(eval(row[r_idx]))
        cedit = row[cedit_idx]
        if cedit == '':
            continue

        dedit = int(float(row[dedit_idx]))
        dmatch = int(float(row[dmatch_idx]))
        nrepeats_idx['dedit'][n_repeats].append(dedit)
        totals['dedit'].append(dedit)

        dq = -10 * np.log10(1 - dmatch / (dmatch + dedit))
        nrepeats_idx['dq'][n_repeats].append(dq)
        totals['dq'].append(dq)

        cedit = int(float(row[cedit_idx]))
        cmatch = int(float(row[cmatch_idx]))

        cq = -10 * np.log10(1 - cmatch / (cmatch + cedit))

        nrepeats_idx['cedit'][n_repeats].append(cedit)
        totals['cedit'].append(cedit)
        nrepeats_idx['cq'][n_repeats].append(cq)
        totals['cq'].append(cq)

    print(np.std(totals['dedit']))
    print(np.std(totals['cedit']))
    plt.hist(totals['dedit'], alpha=0.5, label='SSSL', range=(0, 100), bins=100)
    plt.hist(totals['cedit'], alpha=0.5, label='Consensus', range=(0, 100), bins=100)
    plt.xlim(0, 100)
    plt.legend()
    plt.xlabel("Source Sequence Edit Distance")
    plt.ylabel("# of reads")
    plt.savefig('edit_hist.png')
    plt.clf()

    plt.hist(totals['dq'], alpha=0.5, label='Denoised', range=(0, 40), bins=40)
    plt.hist(totals['cq'], alpha=0.5, label='Consensus', range=(0, 40), bins=40)
    plt.xlim(0, 40)
    plt.legend()
    plt.xlabel("Scaffold Q Concordance")
    plt.ylabel("# of reads")
    plt.savefig('q_hist.png')
    plt.clf()

    cedit_keys = nrepeats_idx['cedit'].keys()
    dedit_keys = nrepeats_idx['dedit'].keys()
    cedit_boxd = [nrepeats_idx['cedit'][k] for k in cedit_keys]
    dedit_boxd = [nrepeats_idx['dedit'][k] for k in dedit_keys]

    fig, ax = plt.subplots()
    bpc = ax.boxplot(cedit_boxd, positions=[r - 0.16 for r in cedit_keys], widths=0.25, patch_artist=True, boxprops=dict(facecolor="C0"), flierprops={'marker': 'o', 'markersize': 3})
    bpd = ax.boxplot(dedit_boxd, positions=[r + 0.16 for r in dedit_keys], widths=0.25, patch_artist=True, boxprops=dict(facecolor="C2"), flierprops={'marker': 'o', 'markersize': 3})
    ax.legend([bpc['boxes'][0], bpd['boxes'][0]], ['Consensus', 'Denoised'], loc='upper right')
    keys = list(range(max(max(cedit_keys), max(dedit_keys))))
    ax.set_xticks(keys, [str(k) for k in keys])
    ax.set_xlabel("# of Subreads")
    ax.set_xlim(1, 20)
    ax.set_ylabel("Source Sequence Edit Distance")

    plt.savefig('edit_box.png')


    print(np.average(lens))
    print(np.average(ns))
    print(min(ns))
    print(np.average(totals['dedit']))
    print(np.average(totals['cedit']))
    cavgs = {}
    cstds = {}
    davgs = {}
    dstds = {}
    for k in nrepeats_idx['dedit']:
        davgs[k] = np.average(nrepeats_idx['dedit'][k])
        dstds[k] = np.std(nrepeats_idx['dedit'][k])
    for k in nrepeats_idx['cedit']:
        cavgs[k] = np.average(nrepeats_idx['cedit'][k])
        cstds[k] = np.std(nrepeats_idx['cedit'][k])

    print('denoising')
    for k in np.sort(list(davgs.keys())):
        print('{}, {}, {}'.format(k, davgs[k], dstds[k]))
    print('consensus')
    for k in np.sort(list(cavgs.keys())):
        print('{}, {}, {}'.format(k, cavgs[k], cstds[k]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='generation file', type=str, required=True)
    parser.add_argument('-n', '--name', help='data name', type=str, required=True)
    args = parser.parse_args()
    gen_plots(args.input, args.name)
