from argparse import ArgumentParser
import os
import pandas as pd

metrics = ['Total Sequences (M)', 'Sequence length (bp)', '%GC', 'Total assembly size (Mbp)', 'Largest contig (Kbp)', 'N50 (Kbp)']


def qc_to_dict(filepath):
    descriptions = {}

    with open(filepath, 'r') as file:
        for line in file.readlines():
            line = line.rstrip()
            metric, expected, value, quality = line.split('\t')
            descriptions[metric] = f"{quality} ({value})"
    
    return descriptions

def read_stadard_qc(dir):
    qc = {}

    for filename in os.listdir(dir):
        if filename.endswith('QC_summary_table.tsv'):
            file_path = os.path.join(dir, filename)

            id = filename.replace('_QC_summary_table.tsv', '')
            qc[id] = qc_to_dict(file_path)
    return pd.DataFrame.from_dict(qc, orient='index').reset_index().rename(columns={'index': 'id'})

def read_bracken(dir):
    dfs = []
    for filename in os.listdir(dir):
        if filename.endswith('_bracken_summary_table.tsv'):
            file_path = os.path.join(dir, filename)

            id = filename.replace('_bracken_summary_table.tsv', '')
            df = pd.read_csv(file_path, sep='\t')
            df['id'] = id
            dfs.append(df)

    df = pd.concat(dfs)
    return df

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-d', '--dir', required=True,
                        help='Directory containing qc tsv')
    parser.add_argument('-o', '--output', required=True,
                        help='Ouput file')
    args = parser.parse_args()
    dir = args.dir
    output_file = args.output

    qc_df = read_stadard_qc(dir)
    bracken_df = read_bracken(dir)

    df = qc_df.merge(bracken_df, on='id', how='left')
    df.to_csv(output_file, index=False)