#! /usr/bin/env python3

from argparse import ArgumentParser, SUPPRESS
from pathlib import Path
import datetime
import logging
import json
import jsonschema
import csv
import re

def process_toxin_coding_genes_args(parser):
    parser.add_argument('-c', '--catalogue', required=True,
                            help='Path to toxin coding gene catalogue')
    parser.add_argument('-s', '--schema', required=False, default=None,
                            help='Path to resistance catalogue')
    parser.add_argument('-b', '--blast_output_tsv', required=True,
                            help='Path to blast output tsv of gene information')
    parser.add_argument('-o', '--output_json', required=False, default="toxin_coding_genes_report.json",
                            help='Path to output json')
    return parser

def initialise_toxin_dict(catalogue:dict):
    toxin_dict = dict()
    for gene in catalogue["toxin_genes"]:
        toxin_dict[gene] = {"presence": False, "percent_identity": -1, "length": -1}
    return toxin_dict

def search_catalogue(catalogue:dict, feature_list:set, toxin_dict:dict):
    for gene in catalogue["toxin_genes"]:
        gene_regex = fr"{gene}"
        for feature_entry in feature_list:
            if re.search(gene_regex, feature_entry["sseqid"]):
                if toxin_dict[gene]["presence"]:
                    if int(feature_entry["length"]) < toxin_dict[gene]["length"]:
                        pass
                    elif int(feature_entry["length"]) == toxin_dict[gene]["length"]:
                        if feature_entry["pident"] > toxin_dict[gene]["pident"]:
                            toxin_dict[gene] = {"presence": True, "percent_identity": float(feature_entry["pident"]), "length": int(feature_entry["length"])}
                    else:
                        toxin_dict[gene] = {"presence": True, "percent_identity": float(feature_entry["pident"]), "length": int(feature_entry["length"])}
                else:
                    toxin_dict[gene] = {"presence": True, "percent_identity": float(feature_entry["pident"]), "length": int(feature_entry["length"])}
    return toxin_dict

def process_toxin_coding_genes(blast_output_tsv: str, catalogue_file: str, schema_file:str, output_json:str):
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

    blast_list = []
    # load sample(s) TSVs
    if not Path(blast_output_tsv).is_file():
        logging.error("{} is not a file.".format(blast_output_tsv))
        raise FileNotFoundError
    else:
        with open(blast_output_tsv) as file:
            reader = csv.DictReader(file, delimiter="\t")
            # build dict of genes/ alleles
            for line in reader:
                blast_list.append(line)
    
    toxin_dict = initialise_toxin_dict(catalogue)
    toxin_dict = search_catalogue(catalogue, blast_list, toxin_dict)

    # Serializing and outputting json
    drug_json = json.dumps(toxin_dict, indent=4)
    with open(output_json, "w") as outfile:
        outfile.write(drug_json)

if __name__ == "__main__":
    parser = ArgumentParser(description='Process BLAST and AMR finder plus outputs to produce AMR report JSON')
    parser = process_toxin_coding_genes_args(parser)
    args = parser.parse_args()

    logging.basicConfig(handlers=[
            #logging.FileHandler('/var/log/GPAS_fetch/{}-{:%Y-%m-%d_%H:%M:%S}-CLIMB-fetch.log'.format(args.run_name, datetime.now())),
            #logging.FileHandler('{}/cdiff_reporting/process_AMR-{:%Y-%m-%d_%H:%M:%S}.log'.format(datetime.now())),
            logging.StreamHandler()],
        format='%(asctime)s - %(levelname)s - %(message)s', 
        level=logging.INFO)
    
    process_toxin_coding_genes(args.blast_output_tsv, args.catalogue, args.schema, args.output_json)