import argparse
import matplotlib.cm as cm
import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 13})
from collections import defaultdict
from scipy.stats import kendalltau, pearsonr, spearmanr, ttest_rel

#slopes = [66.7685891735374, 181.07560740646284, 237.83790970206817, 309.9269443012024, 331.3032814057067, 378.0386443870583, 390.5329308096877, 401.6202871919503, 421.9523129422472, 429.00271269734714, 404.0556795736824, 440.82015852187754, 432.85697478225944, 417.12684942084655, 444.2700164123452, 461.0699080384541, 487.9622702186482]
#slopes = [60.0914587870019, 91.19845874018745, 177.0315827936783, 154.69475235696086, 255.2609291510502, 308.7069323056577, 433.4932813228716, 370.5131280459886, 513.0866466490971, 587.4124314897313, 609.1229993328179, 522.2216469785715, 578.6573301217502, 680.7936516802938, 697.1086431373335, 794.6141519369519, 857.8441198533278]
slopes = [98.77449616346502, 233.90777097785943, 300.8673201068263, 383.74931501139594, 407.9982135943618, 469.1193079395233, 491.0757383674549, 507.6462865595412, 532.3960614615677, 559.0402135317686, 542.135294995765, 599.4681142284905, 569.163216833661, 582.0229732827099, 581.587079112429, 652.9435574357043, 694.099296931412]
#slopes = [242.16190513] * len(slopes)

def gen_figs(fname, psr, k, beta):
    if psr:
        res = json.load(open(fname))
        ncpoints = defaultdict(list)
        ndpoints = defaultdict(list)
        for row in res:
            nr = row[1]
            dent = row[4][str(beta)]
            cent = row[3][str(beta)]
            ncpoints[nr].append(cent)
            ndpoints[nr].append(dent)
        fig = plt.figure()
        ax = fig.add_subplot(111)

        colors = cm.plasma(np.linspace(0, 1, 31))
        diffs = []
        tdiffs = []
        color = []
        totals = 0
        ltotal = 0
        for nr in range(2, 31):
            if nr in ncpoints:
                diffs.append([c - d for c, d in zip(ncpoints[nr], ndpoints[nr])])
                #print(ttest_rel(ndpoints[nr], ncpoints[nr], alternative='less'))
                #print(np.average(diffs[-1]) * slopes[min(nr - 2, len(slopes) - 1)])
                totals += np.sum(diffs[-1]) * slopes[min(nr - 2, len(slopes) - 1)]
                tdiffs.extend(diffs[-1])
                #print(','.join([str(nr), str(len(tdiffs)/len(res))]))
                #print(','.join([str(nr), str(len(diffs))]))
                color.append(colors[nr])
        print(len(tdiffs))

        print(np.average([x for v in ncpoints.values() for x in v]))
        print(np.average([x for v in ndpoints.values() for x in v]))

        total_avg = np.sum(np.sum(diffs)) / np.sum([len(d) for d in diffs])
        print(totals / np.sum([len(d) for d in diffs]))

        print(np.sum(np.sum(diffs))/(np.sum([len(d) for d in diffs])))
        ax.hist(diffs, bins=30, color=color, stacked=True, range=[-0.05, 0.1])
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('plasma'))
        sm.set_clim(vmin=0, vmax=30)
        ax.set_xlim([-0.05, 0.1])
        ax.axvline(total_avg, color='black', linestyle='--', alpha=0.8, label='Average Diff')
        ax.axvline(0, color='black', alpha=0.8)
        plt.legend()
        cbar = plt.colorbar(sm)
        cbar.set_label("# Subreads")
        ax.set_xlabel('Consensus Entropy - SSSL Entropy')
        ax.set_ylabel("# Reads")
        plt.gcf().subplots_adjust(left=0.15)
        fig.savefig('test.png')

        for nr in range(2, 31):
            diffs = [c - d for c, d in zip(ncpoints[nr], ndpoints[nr])]
            diff_mean = np.median(diffs)
            diff_plus = np.quantile(diffs, 0.75)
            diff_minus = np.quantile(diffs, 0.25)
            print(', '.join([str(nr), str(diff_mean), str(diff_plus), str(diff_minus)]))


    else:
        res = json.load(open(fname))
        ncpoints = defaultdict(list)
        ndpoints = defaultdict(list)
        bad = 0
        total = 0
        for row in res:
            nr = row[1]
            if nr == 1:
                continue
            sent = row[2][str(beta)][k-2]
            dent = row[4][str(beta)][k-2]
            if row[3] is not None:
                cent = row[3][str(beta)][k-2]
                ncpoints[nr].append([cent - sent, int(row[-2])])
                if dent < cent and row[-2] < row[-1]:
                    bad += 1
                total += 1
            ndpoints[nr].append([dent - sent , int(row[-1])])
        print(bad/total)
        ws = []
        colors = cm.plasma(np.linspace(1, 0, len(ndpoints.keys())))[::-1]
        for nr, c in zip(sorted(ndpoints.keys()), colors):
            dxs = [p[0] for p in ndpoints[nr]]
            dys = [p[1] for p in ndpoints[nr]]
            plt.scatter(dxs, dys, color=c, alpha=0.5, s=9)
            if nr in ncpoints:
                cxs = [p[0] for p in ncpoints[nr]]
                cys = [p[1] for p in ncpoints[nr]]
                plt.scatter(cxs, cys, color=c, alpha=0.5, s=9)
                dxs.extend(cxs)
                dys.extend(cys)
            print(', '.join([str(nr), str(pearsonr(dxs, dys)[0])]))
            #points = np.linspace(min(dxs), max(dxs), 100)
            #plt.plot(points, np.poly1d(np.polyfit(dxs, dys, 1))(points), color=c)

        plt.xlabel("Source Sequence Entropy Difference")
        plt.ylabel("Source Sequence Edit Distance")
        plt.clim(1, 20)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('plasma'))
        sm.set_clim(vmin=0, vmax=20)
        cbar = plt.colorbar(sm)
        cbar.set_label("# Subreads")
        #plt.plot(np.unique(totaldx + totalcx), np.poly1d(np.polyfit(totaldx + totalcx, totaldy + totalcy, 1))(np.unique(totaldx + totalcx)))
        plt.ylim(bottom=0)
        plt.savefig('test.png')

        plt.clf()
        total_diffs = []
        total_edits = []
        a_s = []
        colors = cm.plasma(np.linspace(0, 1, len(ncpoints.keys())))
        for nr, c in zip(sorted(ncpoints.keys()), colors):
            xs = [c[0] - d[0] for d, c in zip(ndpoints[nr], ncpoints[nr])]
            total_diffs.extend(xs)
            ys = [c[1] - d[1] for d, c in zip(ndpoints[nr], ncpoints[nr])]
            total_edits.extend(ys)
            plt.scatter(xs, ys, color=c, alpha=0.5, s=9)
            xs = np.asarray(xs)
            ys = np.asarray(ys)
            nxs = xs[:, np.newaxis]
            a, _, _, _ = np.linalg.lstsq(nxs, ys)
            a_s.append(a[0])
            yhats = a * xs
            r2 = 1 - np.sum((ys - yhats)**2)/np.sum((ys - np.mean(ys))**2)
            #print(', '.join([str(nr), str(pearsonr(xs, ys)[0])]))
            print(', '.join([str(nr), str(np.sqrt(r2))]))
        print(a_s)
        print(pearsonr(total_diffs, total_edits)[0])
        total_diffs = np.asarray(total_diffs)
        total_edits = np.asarray(total_edits)
        total_diffsn = total_diffs[:, np.newaxis]
        a, _, _, _  = np.linalg.lstsq(total_diffsn, total_edits)
        print(a)
        total_ehats = a * total_diffs
        r2 = 1 - np.sum((total_edits - total_ehats)**2)/np.sum((total_edits - np.mean(total_ehats))**2)
        print(np.sqrt(r2))
        print(np.average(total_diffs))
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap('plasma'))
        sm.set_clim(vmin=0, vmax=20)
        cbar = plt.colorbar(sm)
        cbar.set_label("# Subreads")
        plt.axvline(0, color='black', linestyle='--', alpha=0.3)
        plt.axhline(0, color='black', linestyle='--', alpha=0.3)
        plt.xlabel("Difference in Entropy")
        plt.ylabel("Difference in Source Edit Distance")
        plt.savefig('test1.png')
        #plt.scatter([d - c for d, c in zip(totaldx, totalcx)], [d - c for d, c in zip(totaldy, totalcy)], color='blue', s=9)
        #plt.savefig('test1.png')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='json entropy file')
    parser.add_argument('-k', help='kernel k value', type=int, default=6)
    parser.add_argument('-b', '--beta', help='kernel beta value', type=float, default=8.0)
    parser.add_argument('--psr', action='store_true')
    args = parser.parse_args()
    gen_figs(args.file, args.psr, args.k, args.beta)
