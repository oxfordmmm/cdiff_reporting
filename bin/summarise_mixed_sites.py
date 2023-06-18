from argparse import ArgumentParser
import os
import pandas as pd

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-d', '--dir', required=True,
                        help='Directory containing qc tsv')
    parser.add_argument('-o', '--outfile_prefix', required=True,
                        help='Ouput file')
    args = parser.parse_args()
    dir = args.dir
    outfile_prefix = args.outfile_prefix

    dfs = []
    for filename in os.listdir(dir):
        if filename.endswith('_mixed_sites_summarised.csv'):
            file_path = os.path.join(dir, filename)
            id = filename.replace('_mixed_sites_summarised.csv', '')
            
            df = pd.read_csv(file_path)
            df.insert(0, 'id', id)
            dfs.append(df)
    
    df = pd.concat(dfs)
    df.to_csv(f"{outfile_prefix}_summaries.csv", index=False)
    
    dfs = []
    for filename in os.listdir(dir):
        if filename.endswith('_mixed_sites.csv'):
            file_path = os.path.join(dir, filename)
            id = filename.replace('_mixed_sites.csv', '')
            
            df = pd.read_csv(file_path)
            df.insert(0, 'id', id)
            dfs.append(df)
    
    df = pd.concat(dfs)
    df.to_csv(f"{outfile_prefix}_variant_sites.csv", index=False)