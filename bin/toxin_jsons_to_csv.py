import json
import csv
from argparse import ArgumentParser
import os

def json_to_csv(json_filepath, filename, csv_writer):
    with open(json_filepath, 'r') as json_file:
        data = json.load(json_file)

    # Write data rows
    for gene, gene_data in data.items():
        presence = gene_data['presence']
        percent_identity = gene_data['percent_identity']
        length = gene_data['length']
        gene_length = gene_data['gene_length']
        csv_writer.writerow([filename, gene, presence, percent_identity, length, gene_length])

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-d', '--dir', required=True,
                        help='Directory containing toxin json files')
    parser.add_argument('-o', '--output', required=True,
                        help='Ouput file')
    args = parser.parse_args()
    json_dir = args.dir
    output_file = args.output

    csv_file = open(output_file, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(['Sample', 'Gene', 'Presence', 'Percent Identity', 'Length', 'Gene Length'])

    for filename in os.listdir(json_dir):
        if filename.endswith('toxin_coding_genes_report.json'):
            file_path = os.path.join(json_dir, filename)
            id = filename.replace('_toxin_coding_genes_report.json', '')
            json_to_csv(file_path, id, csv_writer)
    
    csv_file.close()