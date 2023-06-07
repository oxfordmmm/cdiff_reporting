import csv
from argparse import ArgumentParser
import os

metrics = ['Total Sequences (M)', 'Sequence length (bp)', '%GC', 'Total assembly size (Mbp)', 'Largest contig (Kbp)', 'N50 (Kbp)']

def qc_to_csv(filepath, filename, csv_writer):
    descriptions = {}

    with open(filepath, 'r') as file:
        for line in file.readlines():
            line = line.rstrip()
            metric, expected, value, quality = line.split('\t')
            descriptions[metric] = f"{quality} ({value})"
    
    csv_writer.writerow([
        id,
        descriptions[metrics[0]],
        descriptions[metrics[1]],
        descriptions[metrics[2]],
        descriptions[metrics[3]],
        descriptions[metrics[4]],
        descriptions[metrics[5]]
    ])

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-d', '--dir', required=True,
                        help='Directory containing qc tsv')
    parser.add_argument('-o', '--output', required=True,
                        help='Ouput file')
    args = parser.parse_args()
    dir = args.dir
    output_file = args.output

    csv_file = open(output_file, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(['Sample'] + metrics)

    for filename in os.listdir(dir):
        if filename.endswith('QC_summary_table.tsv'):
            file_path = os.path.join(dir, filename)
            id = filename.replace('_QC_summary_table.tsv', '')
            qc_to_csv(file_path, id, csv_writer)
    
    csv_file.close()