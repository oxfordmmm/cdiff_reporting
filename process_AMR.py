#! /usr/bin/env python3

from argparse import ArgumentParser, SUPPRESS
from pathlib import Path
import datetime
import logging
import json
import jsonschema
import csv
import re

import utils

def process_AMR_args(parser):
    parser.add_argument('-c', '--catalogue', required=True,
                            help='Path to resistance catalogue')
    parser.add_argument('-s', '--schema', required=False, default=None,
                            help='Path to resistance catalogue')
    parser.add_argument('-b', '--blast_output_tsv', required=True,
                            help='Path to blast output tsv of gene information')
    parser.add_argument('-f', '--amr_finder_output_tsv', required=True,
                            help='Path to amr finder plus output tsv of mutation information')
    parser.add_argument('-o', '--output_json', required=False, default="resistance_report.json",
                            help='Path to output json')
    return parser


def process_AMR(blast_output_tsv: str, amr_finder_output_tsv:str, catalogue_file: str, schema_file:str, output_json:str):
    # load JSON catalogue
    if not Path(catalogue_file).is_file():
        logging.error("{} is not a file.".format(catalogue_file))
        raise FileNotFoundError
    
    if schema_file and not Path(schema_file).is_file():
        logging.error("{} is not a file.".format(schema_file))
        raise FileNotFoundError

    with open(catalogue_file) as f:
        catalogue = json.loads(f.read())
    
    with open(schema_file) as f:
        schema = json.loads(f.read())

    jsonschema.validate(instance=catalogue, schema=schema)

    blast_list = set()
    amr_finder_list = set()
    # load sample(s) TSVs
    if not Path(blast_output_tsv).is_file():
        logging.error("{} is not a file.".format(blast_output_tsv))
        raise FileNotFoundError
    else:
        with open(blast_output_tsv) as file:
            reader = csv.reader(file, delimiter="\t")
            # build dict of genes/ alleles
            for line in reader:
                blast_list.add(line[1])
    
    if not Path(amr_finder_output_tsv).is_file():
        logging.error("{} is not a file.".format(amr_finder_output_tsv))
        raise FileNotFoundError
    else:
        with open(amr_finder_output_tsv) as file:
            reader = csv.reader(file, delimiter="\t")
            # build set of mutations
            for line in reader:
                amr_finder_list.add(line[5])
    
    drug_resistances = dict()
    for drug in catalogue["drugs"]:
        drug_resistances[drug] = "S"

    for gene, gene_value in catalogue["genes"].items():
        gene_regex = fr"{gene}"
        for feature in blast_list:
            if re.search(gene_regex, feature):
                if "alleles" in gene_value:
                    for allele, allele_value in gene_value["alleles"].items():
                        allele_regex = fr"{allele}"
                        if re.search(allele_regex, feature):
                            if "mutations" in allele_value:
                                for mutation, mutation_value in allele_value["mutations"].items():
                                    mutation_regex = fr"{mutation}"
                                    if re.search(mutation_regex, feature):
                                        drug_resistances[mutation_value["drug"]] = "R"
                            else:
                                drug_resistances[allele_value["drug"]] = "R"
                else:
                    drug_resistances[gene_value["drug"]] = "R"

    # Serializing and outputting json
    drug_json = json.dumps(drug_resistances, indent=4)
    with open(output_json, "w") as outfile:
        outfile.write(drug_json)

if __name__ == "__main__":
    parser = ArgumentParser(description='Process BLAST and AMR finder plus outputs to produce AMR report JSON')
    parser = process_AMR_args(parser)
    args = parser.parse_args()

    logging.basicConfig(handlers=[
            #logging.FileHandler('/var/log/GPAS_fetch/{}-{:%Y-%m-%d_%H:%M:%S}-CLIMB-fetch.log'.format(args.run_name, datetime.now())),
            #logging.FileHandler('{}/cdiff_reporting/process_AMR-{:%Y-%m-%d_%H:%M:%S}.log'.format(datetime.now())),
            logging.StreamHandler()],
        format='%(asctime)s - %(levelname)s - %(message)s', 
        level=logging.INFO)
    
    process_AMR(args.blast_output_tsv, args.amr_finder_output_tsv, args.catalogue, args.schema, args.output_json)