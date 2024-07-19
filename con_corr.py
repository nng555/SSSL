import csv
from collections import defaultdict
import os
import ast
from scipy.stats import pearsonr
from nltk.metrics import edit_distance
import argparse

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('csv_file', help='The name of the CSV file')
args = parser.parse_args()

# Open the CSV file
with open(args.csv_file, newline='') as csvfile, \
        open(args.csv_file[:-4] + '_new.csv', 'w') as outfile:
    reader = csv.DictReader(csvfile)
    writer = csv.DictWriter(outfile, reader.fieldnames + ['consensus_r_edit', 'denoised_r_edit'])
    writer.writeheader()

    # Initialize variables
    # Loop through each row in the CSV file
    for i, row in enumerate(reader):
        if (i + 1) % 100 == 0:
            print(i + 1)

        noisy_seqs = ast.literal_eval(row['read'])
        nrepeats = len(noisy_seqs)

        row['consensus_r_edit'] = ''
        row['denoised_r_edit'] = ''
        if row['consensus'] != '':

            # Calculate the average edit distance between the consensus and each sequence in noisy_repeats
            noisy_ed = 0
            for seq in noisy_seqs:
                noisy_ed += edit_distance(row['consensus'], seq)
            noisy_ed /= len(noisy_seqs)
            row['consensus_r_edit'] = noisy_ed

        if row['denoised'] != '':

            # Calculate the average edit distance between the consensus and each sequence in noisy_repeats
            noisy_ed = 0
            for seq in noisy_seqs:
                noisy_ed += edit_distance(row['denoised'], seq)
            noisy_ed /= len(noisy_seqs)
            row['denoised_r_edit'] = noisy_ed

        writer.writerow(row)

