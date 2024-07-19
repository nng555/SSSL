import csv
from collections import defaultdict
import os
import ast
from scipy.stats import pearsonr, linregress
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('csv_file', help='The name of the CSV file')
parser.add_argument('base_file', help='name of baseline csv file')
args = parser.parse_args()

base_data = pd.read_csv(args.base_file)
gen_data = pd.read_csv(args.csv_file)
gen_data = gen_data.drop('scaffold', axis=1)

new_data = pd.merge(base_data, gen_data, right_on='read', left_on='noisy_repeats', how='left')
new_fname = args.csv_file.split('.csv')[0] + '_baseline.csv'
new_data.to_csv(new_fname)
