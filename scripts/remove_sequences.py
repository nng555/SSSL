import pandas as pd
import argparse

COLS = [
    'read',
    'noisy_read',
    'denoised',
    'scaffold',
    'consensus',
    'muscle',
    'mafft',
    'tcoffee',
    'median',
    'quickmedian'
]

def remove_sequences(fname):
    dframe = pd.read_csv(fname)
    for col in COLS:
        if col in dframe.columns:
            dframe = dframe.drop(columns=col)

    dframe.to_csv(fname.split('.csv')[0] + '_private.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input file to remove sensitive data')
    args = parser.parse_args()
    remove_sequences(args.input)
