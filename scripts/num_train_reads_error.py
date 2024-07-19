import argparse
import os
import csv

def analyze(train_dir, output):

    train_count = {}
    for f in os.listdir(train_dir):
        if 'data_shard' not in f:
            continue
        with open(os.path.join(train_dir, f)):
            reader = csv.DictReader(

    with open(output, 'r') as


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_arg('-t', '--train-dir', help='directory of training shards')
    parser.add_arg('-o', '--output', help='output generated test file')
    args = parser.parse_args()

    analyze(args.train_dir, args.output)
