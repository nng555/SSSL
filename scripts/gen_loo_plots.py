import csv
from matplotlib import ticker
import matplotlib.cm as cm
from collections import defaultdict
import ast
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
from scipy.stats import pearsonr, ttest_rel

def gen_loo_plots(cfname, gfname):

    def load_edits(fname, edit_name, res):
        print(edit_name)
        print(fname)
        reader = csv.reader(open(fname))
        header = next(reader)
        r_idx = header.index('read')
        d_idx = header.index('denoised')
        ce_idx = header.index('consensus_r_edit')
        de_idx = header.index('denoised_r_edit')
        l_idx = header.index('loo_edit')
        for row in reader:
            skey = row[r_idx]
            if row[l_idx] == '':
                loo_edit = None
            elif '[' in row[l_idx]:
                tmp = ast.literal_eval(row[l_idx])
                if any([isinstance(l, list) for l in tmp]):
                    loo_edit = np.average(tmp[0])
                else:
                    loo_edit = np.average(tmp)
            else:
                loo_edit = int(row[l_idx])

            if skey in res:
                if edit_name + '_loo_edit' in res[skey]:
                    #raise Exception("Key clash")
                    continue
                #assert res[skey]['denoised_edit'] == row[de_idx]
                #assert res[skey]['consensus_edit'] == row[ce_idx]
                res[skey][edit_name + '_loo_edit'] = loo_edit
            else:
                cedit = int(row[ce_idx]) if row[ce_idx] != "" else None

                res[skey] = {
                    'consensus_edit': cedit,
                    'denoised_edit': int(row[de_idx]),
                    edit_name + '_loo_edit': loo_edit,
                    'n_repeats': row[r_idx].count(',') + 1,
                }
        return res

    res = load_edits(cfname, 'consensus', {})
    res = load_edits(gfname, 'denoised', res)

    dframe = defaultdict(list)
    for v in res.values():
        for k in v.keys():
            dframe[k].append(v[k])

    dframe = pd.DataFrame(dframe)
    dframe['edit_diffs'] = dframe['consensus_edit'] - dframe['denoised_edit']
    dframe['loo_diffs'] = dframe['consensus_loo_edit'] - dframe['denoised_loo_edit']

    dmask = dframe['denoised_edit'].notna() & dframe['denoised_loo_edit'].notna()
    cmask = dframe['consensus_edit'].notna() & dframe['consensus_loo_edit'].notna()
    tmask = dmask & cmask

    print(pearsonr(dframe['edit_diffs'][tmask], dframe['loo_diffs'][tmask]))
    plt.scatter(dframe['edit_diffs'][tmask].to_numpy(), dframe['loo_diffs'][tmask].to_numpy())
    plt.show()
    plt.savefig('diff_scatter.png')
    plt.cla()

    plt.scatter(dframe['consensus_edit'][cmask], dframe['consensus_loo_edit'][cmask])
    plt.scatter(dframe['denoised_edit'][dmask], dframe['denoised_loo_edit'][dmask])
    plt.show()
    plt.savefig("sep_scatter.png")
    plt.cla()

    keys = np.sort(dframe['n_repeats'][tmask].unique())

    colors = cm.plasma(np.linspace(0, 1, len(keys)))
    for nr, c in zip(keys, colors):
        kmask = dframe['n_repeats'] == nr
        plt.scatter(dframe['loo_diffs'][tmask & kmask], dframe['edit_diffs'][tmask & kmask], color=c, alpha=0.5, s=9)
        #print(','.join([str(v) for v in [nr, np.average(dframe['loo_diffs'][tmask & kmask]), np.std(dframe['loo_diffs'][tmask & kmask])]]))
        print(','.join([str(v) for v in [nr, np.median(dframe['loo_diffs'][tmask & kmask]), np.quantile(dframe['loo_diffs'][tmask & kmask], 0.75), np.quantile(dframe['loo_diffs'][tmask & kmask], 0.25)]]))
        print(len(dframe['loo_diffs'][tmask & kmask]))
        #print(ttest_rel(dframe['consensus_loo_edit'][tmask & kmask], dframe['denoised_loo_edit'][tmask & kmask], alternative='greater'))
        #print((np.average(dframe['loo_diffs'][tmask & kmask]), np.average(dframe['edit_diffs'][tmask & kmask])))
        #print(np.average((dframe['loo_diffs'] - dframe['edit_diffs'])[tmask & kmask]))

    xpoints = ypoints = plt.xlim()
    #plt.plot(xpoints, ypoints,  color='k', scalex=False, scaley=False, alpha=0.6)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('plasma'))
    sm.set_clim(vmin=2, vmax=20)
    cbar = plt.colorbar(sm)
    cbar.ax.set_yticklabels([v for v in range(2, 21, 2)])
    cbar.set_label("# Subreads")
    plt.axvline(0, color='black', linestyle='--', alpha=0.3)
    plt.axhline(0, color='black', linestyle='--', alpha=0.3)
    plt.xlabel("Difference in LOO Edit")
    plt.ylabel("Difference in Source Edit Distance")
    plt.savefig('test1.png')

    for k in range(3, 20):
        kmask = dframe['n_repeats'] == k
        print(','.join([str(k), str(pearsonr(dframe['edit_diffs'][tmask & kmask], dframe['loo_diffs'][tmask & kmask])[0])]))

    print(np.average(dframe['consensus_loo_edit'][tmask]))
    print(np.average(dframe['denoised_loo_edit'][tmask]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--consensus', help='consensus loo file path')
    parser.add_argument('-g', '--generate', help='generated sssl loo file path')
    args = parser.parse_args()
    gen_loo_plots(args.consensus, args.generate)
