#! /env/bin/python3
 
from argparse import ArgumentParser, SUPPRESS
import logging
from pathlib import Path
import json
import csv

def make_summary_json_args(parser):
    parser.add_argument('-i', '--sample_list', required=True,
                            help='Path to sample list used for runlistcompare')
    parser.add_argument('-c', '--cluster_dir', required=True,
                            help='Path to directory for cluster.txt and trees dir')
    parser.add_argument('-d', '--dirs',  nargs="*", required=True,
                            help='Path to directories that will be searched for QC files')
    parser.add_argument('-o', '--output_json', required=False, default="summary_data.json",
                            help='Path to output JSON')
    return parser

def make_summary_json(sample_list: str, dirs:list[str], cluster_dir:str, output_json:str):
    samples_json = {}
    with open(sample_list) as file:
        for line in csv.reader(file, delimiter="\t"):
            samples_json[line[0]] = {"qc_file": ""}
            for dir in dirs:
                for match in Path(dir).glob(f'{line[0]}*QC_summary_table.tsv'):
                    samples_json[line[0]]["qc_file"] = match
                    break
    
    cluster_path = Path(cluster_dir)
    samples_json["cluster_file"] = str(cluster_path / "clusters.txt")
    samples_json["cluster_trees_dir"] = str(cluster_path / "cluster_ml")
    
    with open(output_json, "w") as output_file:
        json.dump(samples_json, output_file, indent=4, default=str)


if __name__ == '__main__':
    parser = ArgumentParser(description='Make the JSON to use for generating a summary report.')
    parser = make_summary_json_args(parser)
    args = parser.parse_args()

    logging.basicConfig(handlers=[
            logging.StreamHandler()],
        format='%(asctime)s - %(levelname)s - %(message)s', 
        level=logging.INFO)

    make_summary_json(args.sample_list, args.dirs, args.cluster_dir, args.output_json)

