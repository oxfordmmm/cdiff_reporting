#! /env/bin/python3
from argparse import ArgumentParser
from os.path import exists
import pandas as pd

def get_bracken_stats(bracken_file, qc_dict):
    print(bracken_file)
    df = pd.read_csv(bracken_file, sep='\t').sort_values('fraction_total_reads', ascending=False, ignore_index=True)
    df = df.iloc[:3]
    top_hits = [(name, round(100*pc,2)) for name, pc in zip(df['name'], df['fraction_total_reads'])]
    qc_dict['bracken_top_3'] = top_hits


def write_qc(qc_dict, outfile):
    with open(outfile, 'w') as file:
        file.write("First\tSecond\tThird\n")
        bracken_strs = [f"{species} ({pc})" for species, pc in qc_dict['bracken_top_3']]
        file.write('\t'.join(bracken_strs) + '\n')



if __name__ == '__main__':
    parser = ArgumentParser(description='Summarise bracken outputs to tsv.')
    parser.add_argument('-b', '--bracken_output', required=True)
    parser.add_argument('-o', '--outfile', required=True)
    args = parser.parse_args()

    qc_dict = {
        'bracken_top_3': [('', 0), ('', 0), ('', 0)]
    }

    if exists(args.bracken_output):
        get_bracken_stats(args.bracken_output, qc_dict)
    else:
        print("bracken file doesn't exist!")

    write_qc(qc_dict, args.outfile)

