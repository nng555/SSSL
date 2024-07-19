from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo
from nltk.metrics import edit_distance

from kmer_trie import KmerTrie

import ast
import numpy as np
import subprocess
import argparse
import csv
import os
import tempfile
import Levenshtein
import itertools

MUSCLE_PATH = '/h/nng/programs/msa/muscle'
TCOFFEE_PATH = '/h/nng/programs/msa/t_coffee'
MAFFT_PATH = '/h/nng/programs/msa/mafft/mafft.bat'

def run_baselines(in_file, num_shards, shard, private):

    def align(bname, infile):
        outfile = tempfile.NamedTemporaryFile()
        if bname == 'muscle':
            p = subprocess.run([MUSCLE_PATH, "-align", infile, '-output', outfile.name], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        elif bname == 'tcoffee':
            p = subprocess.run([TCOFFEE_PATH, '-output', 'fasta', '-infile', infile, '-outfile', outfile.name], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        elif bname == 'mafft':
            p = subprocess.run([MAFFT_PATH, '--localpair', '--maxiterate', '100', infile], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
            with open(outfile.name, 'wb') as mafftf:
                mafftf.write(p.stdout)
        align = AlignIO.read(outfile.name, 'fasta')
        msa_consensus = str(AlignInfo.SummaryInfo(align).gap_consensus(threshold=0)).upper().replace('-', '')
        return msa_consensus

    def steiner(bname, seqs):
        if bname == 'median':
            return Levenshtein.median(seqs)
        elif bname == 'quickmedian':
            return Levenshtein.quickmedian(seqs)

    with open(in_file) as inf, \
         open(in_file.split('.csv')[0] + '_baseline_{}.csv'.format(str(shard)), 'w') as outf:

        print("Begin processing shard {}/{}".format(str(shard + 1), str(num_shards)), flush=True)
        reader = csv.DictReader(inf)

        private_headers = ['noisy_repeats', 'scaffold']
        baselines = ['muscle', 'tcoffee', 'mafft']
        steiners = ['median', 'quickmedian']
        fields = ['read_id', 'nrepeats', 'scaffold_subread_edit']
        for b in baselines + steiners:
            if not private:
                fields.append(b)
            fields += [b + '_edit']
            fields += [b + '_subread_edit', b + '_loo_edit']

        fieldnames = reader.fieldnames + fields
        if private:
            fieldnames = [field for field in fieldnames if field not in private_headers]

        writer = csv.DictWriter(outf, fieldnames=fieldnames)
        writer.writeheader()

        for i, row in enumerate(reader):
            if i % num_shards != shard:
                continue

            if (i // num_shards + 1) % 10 == 0:
                print("Processing row {}, chunk {}".format(str(i + 1), str(i // num_shards + 1)), flush=True)

            noisy_seqs = ast.literal_eval(row['noisy_repeats'])
            row['nrepeats'] = len(noisy_seqs)

            if row['nrepeats'] < 3:
                writer.writerow(row)
                continue

            scaffold = row['scaffold'].upper()
            tmp_if_main = tempfile.NamedTemporaryFile()
            tmp_if_loo = tempfile.NamedTemporaryFile()
            #etrie = KmerTrie(8, unk_char='X')

            with open(tmp_if_main.name, 'w') as f:
                for snum, s in enumerate(noisy_seqs):
                    f.write('>s{}\n{}\n'.format(str(snum), s))
                    #etrie.insert(s, value= 1 / len(noisy_seqs))

            #etrie.insert(scaffold)
            #scaffold_entropy = etrie.get_entropy([8])[8.0]
            #etrie.delete(scaffold)
            #row['scaffold_entropy'] = scaffold_entropy

            scaffold_subread_edit = 0
            for s in noisy_seqs:
                scaffold_subread_edit += edit_distance(scaffold, s)
            scaffold_subread_edit /= len(noisy_seqs)
            row['scaffold_subread_edit'] = scaffold_subread_edit

            for b in steiners + baselines:

                if b in baselines:
                    msa_consensus = align(b, tmp_if_main.name)
                else:
                    msa_consensus = steiner(b, noisy_seqs)
                msa_edit = edit_distance(msa_consensus, scaffold)
                row[b + '_edit'] = msa_edit
                if not private:
                    row[b] = msa_consensus

                # LOO alignments
                loo_edit = 0
                for lidx in range(len(noisy_seqs)):
                    if b in baselines:
                        with open(tmp_if_loo.name, 'w') as f:
                            for j, s in enumerate(noisy_seqs):
                                if lidx == j:
                                    continue
                                f.write('>s{}\n{}\n'.format(str(j), s))
                        loo_consensus = align(b, tmp_if_loo.name)
                    else:
                        loo_consensus = steiner(b, noisy_seqs[:lidx] + noisy_seqs[lidx+1:])
                    loo_edit += edit_distance(loo_consensus, noisy_seqs[lidx])
                loo_edit /= len(noisy_seqs)

                '''
                # LBO alignments
                nseqs = len(noisy_seqs)
                lbo_edits = [[] for _ in range(nseqs)]
                seq_set = set(np.arange(nseqs))

                for _ in range(100):
                    combo = np.random.choice(nseqs, nseqs, replace=True)
                    test_set = seq_set.difference(set(combo))

                    if b in baselines:
                        with open(tmp_if_loo.name, 'w') as f:
                            for j, sidx in enumerate(combo):
                                f.write('>s{}\n{}\n'.format(str(j), noisy_seqs[sidx]))
                        lbo_consensus = align(b, tmp_if_loo.name)
                    else:
                        lbo_consensus = steiner(b, [noisy_seqs[j] for j in combo])
                    for tidx in test_set:
                        lbo_edits[tidx].append(edit_distance(lbo_consensus, noisy_seqs[tidx]))

                lbo_edit = np.mean([np.mean(lbo) for lbo in lbo_edits])

                # LHO alignments
                lho_edit = 0
                nseqs = len(noisy_seqs)
                nchoices = int(np.ceil(nseqs/2))
                if len(noisy_seqs) >= 7:
                    combos = [np.random.choice(np.arange(nseqs), nchoices, replace=False) for _ in range(30)]
                else:
                    combos = list(itertools.combinations(np.arange(nseqs), nchoices))

                nseqs = len(noisy_seqs)
                seq_set = set(np.arange(nseqs))
                for combo in combos:
                    test_set = seq_set.difference(set(combo))
                    if b in baselines:
                        with open(tmp_if_loo.name, 'w') as f:
                            for j in combo:
                                f.write('>s{}\n{}\n'.format(str(j), noisy_seqs[j]))
                        lho_consensus = align(b, tmp_if_loo.name)
                    else:
                        lho_consensus = steiner(b, [noisy_seqs[j] for j in combo])
                    combo_edit = 0
                    for tidx in test_set:
                        combo_edit += edit_distance(lho_consensus, noisy_seqs[tidx])
                    combo_edit /= len(test_set)
                    lho_edit += combo_edit
                lho_edit /= len(combos)
                '''

                #etrie.insert(msa_consensus)
                #entropy = etrie.get_entropy([8])[8.0]
                #etrie.delete(msa_consensus)

                msa_subread_edit = 0
                for s in noisy_seqs:
                    msa_subread_edit += edit_distance(msa_consensus, s)
                msa_subread_edit /= len(noisy_seqs)
                row.update({
                    b + '_subread_edit': msa_subread_edit,
                    b + '_loo_edit': loo_edit,
                    #b + '_lbo_edit': lbo_edit,
                    #b + '_entropy': entropy,
                })

            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input csv file to process')
    parser.add_argument('-n', '--num-shards', help='total number of shards', type=int)
    parser.add_argument('-s', '--shard', help='shard number', type=int)
    parser.add_argument('-p', '--private', help='dont dump sequences', action='store_true')
    args = parser.parse_args()
    run_baselines(args.input, args.num_shards, args.shard, args.private)
