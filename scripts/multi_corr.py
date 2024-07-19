import csv
from collections import defaultdict
import os
import ast
from scipy.stats import pearsonr
from nltk.metrics import edit_distance
import argparse
from multiprocessing import Process, Manager

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('csv_dir', help='The name of the CSV directory')
args = parser.parse_args()

def calculate_edit_distance(csv_file, edit_dist_scaffold, edit_dist_noisy):
    # Open the CSV file
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        # Initialize variables
        # Loop through each row in the CSV file
        for i, row in enumerate(reader):
            if (i + 1) % 100 == 0:
                print(i + 1)

            if row['consensus'] == '':
                continue
            noisy_seqs = ast.literal_eval(row['read'])
            nrepeats = len(noisy_seqs)

            # Calculate the edit distance from the scaffold to the consensus
            scaffold_ed = edit_distance(row['scaffold'], row['consensus'])
            edit_dist_scaffold[nrepeats].append(scaffold_ed)

            # Calculate the average edit distance between the consensus and each sequence in noisy_repeats
            noisy_ed = 0
            for seq in noisy_seqs:
                noisy_ed += edit_distance(row['consensus'], seq)
            noisy_ed /= len(noisy_seqs)
            edit_dist_noisy[nrepeats].append(noisy_ed)

if __name__ == '__main__':
    # Use a manager to share the dictionaries between processes
    with Manager() as manager:
        edit_dist_scaffold = manager.dict()
        edit_dist_noisy = manager.dict()

        # Create a separate process for each CSV file
        processes = []
        for f in os.listdir(args.csv_dir):
            if 'data_shard' not in f:
                continue

            csv_file = os.path.join(args.csv_dir, f)
            p = Process(target=calculate_edit_distance, args=(csv_file, edit_dist_scaffold, edit_dist_noisy))
            processes.append(p)
            p.start()

        # Wait for all processes to finish
        for p in processes:
            p.join()

        # Calculate the Pearson correlation
        for k in edit_dist_scaffold.keys():
            r, p = pearsonr(edit_dist_scaffold[k], edit_dist_noisy[k])

            # Print the result
            print('Pearson correlation for size {}: {}'.format(str(k), str(r)))
