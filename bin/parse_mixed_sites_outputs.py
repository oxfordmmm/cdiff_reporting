#! /env/bin/python3
from argparse import ArgumentParser
import os
import sys
import pandas as pd



if __name__ == '__main__':
    parser = ArgumentParser(description='Parse mixed infection estimate')
    parser.add_argument('-m', '--mixed_infection_estimate', required=True)
    parser.add_argument('-o', '--outfile', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.mixed_infection_estimate):
        print("Error: mixed infection file doesn't exist.")
        sys.exit()

    try:
        df = pd.read_csv(args.mixed_infection_estimate, sep='\t', index_col=False)
        df['mixed_ST'] = df['ML_sites_diff'].iloc[0] > 1.5 and df['deviance'].iloc[0] > 112

        df.to_csv(args.outfile, index=False, sep='\t')
    except Exception as e:
        print("Error occurred while processing mixed infection output:", str(e))


