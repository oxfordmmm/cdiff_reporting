#! /env/bin/python3
from argparse import ArgumentParser
from os.path import exists
import pandas as pd

def get_fastqc_stats(fastqc_file, qc_dict):
    at_length_section = False
    length_sum = 0
    count_sum = 0
    with open(fastqc_file, 'r') as file:
        for line in file.readlines():
            line = line.rstrip()

            if at_length_section:
                if line.startswith('>>END_MODULE'):
                    at_length_section = False
                    continue
                length, count = line.split('\t')
                count = round(float(count))
                min, max = length.split('-')
                length = (int(min) + int(max)) * 0.5
                count_sum += count
                length_sum += count * length
            elif line.startswith('Total Sequences'):
                qc_dict['n_sequences'] = int(line.split('\t')[1])
            elif line.startswith("%GC"):
                qc_dict['read_GC'] = int(line.split('\t')[1])
            elif line.startswith("Sequence length"):
                length = line.split('\t')[1]
                if '-' in length:
                    min, max = length.split('-')
                    qc_dict['seq_length'] = round((int(min) + int(max)) * 0.5)
                else:
                    qc_dict['seq_length'] = int(length)
            elif line.startswith("#Length"):
                at_length_section = True
    
    qc_dict['avg_seq_length'] = round(length_sum/count_sum)


def get_quast_stats(quast_file, qc_dict):
    with open(quast_file, 'r') as file:     
        for line in file.readlines():
            if line.startswith('Total length\t'):
                qc_dict['assembly_size'] = int(line.split('\t')[1])
            elif line.startswith('Largest contig'):
                qc_dict['largest_contig'] = int(line.split('\t')[1])
            elif line.startswith('GC'):
                qc_dict['assembly_GC'] = float(line.split('\t')[1])
            elif line.startswith('N50'):
                qc_dict['N50'] = int(line.split('\t')[1])

def get_bracken_result(bracken_summary, qc_dict):
    df = pd.read_csv(bracken_summary, sep='\t')
    qc_dict['bracken_other_pc'] = df['other_pc'].iloc[0]

def get_mixed_infection_result(mixed_infection_estimate, qc_dict):
    df = pd.read_csv(mixed_infection_estimate, sep='\t', index_col=False)
    qc_dict['mixed_infection'] = df['mixed_ST'].iloc[0]

def get_depth_result(depth, qc_dict):
    df = pd.read_csv(depth)
    qc_dict['avg_depth'] = round(df['avDepth'].iloc[0], 1)

def get_pc_called(n_rate_file, qc_dict):
    results = pd.read_csv(n_rate_file, sep='\t').iloc[0].to_dict()
    qc_dict['pc_called'] = round(results['pc_called'], 2)

def write_qc(qc_dict, outfile):
    with open(outfile, 'w') as file:
        file.write("Metric\tExpected value\tActual value\tQuality\n")

        n_seq = round(qc_dict['n_sequences'] / 1000000, 1)
        quality = 'Pass' if 0.5 <= n_seq <= 20 else 'Fail'
        file.write(f"Total Sequences (M)\t0.5-20\t{n_seq}\t{quality}\n")

        seq_length = qc_dict['avg_seq_length']
        quality = 'Pass' if 50 <= seq_length <= 150 else 'Fail'
        file.write(f"Sequence length (bp)\t50-150\t{seq_length}\t{quality}\n")

        gc = round(qc_dict['assembly_GC'], 1)
        quality = 'Pass' if 27.9 <= gc <= 29.2 else 'Fail'
        file.write(f"%GC\t27.9-29.2\t{gc}\t{quality}\n")

        size = round(qc_dict['assembly_size'] / 1000000, 1)
        quality = 'Pass' if 3.8 <= size <= 4.8 else 'Fail'
        file.write(f"Total assembly size (Mbp)\t3.8-4.8\t{size}\t{quality}\n")

        largest_contig = round(qc_dict['largest_contig'] / 1000, 1)
        quality = 'Pass' if 1 <= largest_contig <= 4300 else 'Fail'
        file.write(f"Largest contig (Kbp)\t1-4300\t{largest_contig}\t{quality}\n")

        n50 = round(qc_dict['N50'] / 1000, 1)
        quality = 'Pass' if 10 <= n50 <= 1000 else 'Fail'
        file.write(f"N50 (Kbp)\t10-1000\t{n50}\t{quality}\n")

        bracken_other = round(qc_dict['bracken_other_pc'], 1)
        quality = 'Pass' if bracken_other <= 2 else 'Fail'
        file.write(f"Residual contamination %\t0-2\t{bracken_other}\t{quality}\n")

        mixed_infection = qc_dict['mixed_infection']
        quality = 'Pass' if not mixed_infection else 'Fail'
        file.write(f"Mixed ST's\tFalse\t{mixed_infection}\t{quality}\n")

        value = qc_dict['avg_depth']
        quality = 'Pass' if value >= 20 else 'Fail'
        file.write(f"avg read depth\t>=20\t{value}\t{quality}\n")

        value = qc_dict['pc_called']
        quality = 'Pass' if value > 95 else 'Fail'
        file.write(f"% of bases called\t>95\t{value}\t{quality}\n")


if __name__ == '__main__':
    parser = ArgumentParser(description='Summarise fastqc and quast outputs to tsv.')
    parser.add_argument('-f', '--fastQC_output', required=True)
    parser.add_argument('-q', '--quast_output', required=True)
    parser.add_argument('-b', '--bracken_summary', required=True)
    parser.add_argument('-m', '--mixed_infection_estimate', required=True)
    parser.add_argument('-d', '--read_depth_stats', required=True)
    parser.add_argument('-n', '--n_count', required=True)
    parser.add_argument('-o', '--outfile', required=True)
    args = parser.parse_args()

    qc_dict = {
        'n_sequences' : 0,
        'seq_length' : 0,
        'avg_seq_length' : 0,
        'read_GC' : 0,
        'assembly_GC' : 0,
        'assembly_size' : 0,
        'largest_contig' : 0,
        'N50' : 0,
        'bracken_other_pc': 0,
        'mixed_infection': False,
        'pc_called' : 0,
        'avg_depth': 0,
    }

    files = [args.fastQC_output, args.quast_output, args.bracken_summary,
             args.mixed_infection_estimate, args.read_depth_stats,
             args.n_count]
    funcs = [get_fastqc_stats, get_quast_stats, get_bracken_result,
             get_mixed_infection_result, get_depth_result,
             get_pc_called]

    for file, func in zip(files, funcs):
        if exists(file):
            func(file, qc_dict)
        else:
            print(f"{file} file doesn't exist, skipping")

    write_qc(qc_dict, args.outfile)

