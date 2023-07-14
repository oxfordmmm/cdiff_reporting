#! /usr/bin/env python3

from argparse import ArgumentParser, SUPPRESS
from pathlib import Path
import datetime
import logging
import json
import jsonschema
import csv
import re
import pandas as pd

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
        gene_length = catalogue["toxin_gene_lengths"][gene]
        toxin_dict[gene] = {"presence": False, "percent_identity": -1,
                            "length": -1, "gene_length": gene_length}
    return toxin_dict

def find_toxin_genes_best(catalogue:dict, feature_list:set, toxin_dict:dict):
    for gene in catalogue["toxin_genes"]:
        gene_regex = fr"{gene}"
        for feature_entry in feature_list:
            if gene_regex == feature_entry["sseqid"]:
                gene_length = catalogue["toxin_gene_lengths"][gene]
                match_pc = 100 * int(feature_entry["length"]) / gene_length
                pident = float(feature_entry["pident"])
                if pident < 90:
                    continue
                if pident < 97 and match_pc < 90:
                    continue
                if toxin_dict[gene]["presence"]:
                    if int(feature_entry["length"]) < toxin_dict[gene]["length"]:
                        pass
                    elif int(feature_entry["length"]) == toxin_dict[gene]["length"]:
                        if pident > toxin_dict[gene]["percent_identity"]:
                            toxin_dict[gene] = {"presence": True, "percent_identity": pident, "length": int(feature_entry["length"]), "gene_length": gene_length}
                    else:
                        toxin_dict[gene] = {"presence": True, "percent_identity": pident, "length": int(feature_entry["length"]), "gene_length": gene_length}
                else:
                    toxin_dict[gene] = {"presence": True, "percent_identity": pident, "length": int(feature_entry["length"]), "gene_length": gene_length}
    return toxin_dict

def find_toxin_genes_sum(catalogue:dict, blast_file:str, toxin_dict:dict):
    toxin_genes = catalogue["toxin_genes"]
    df = pd.read_csv(blast_file, sep='\t')
    df = df.query('pident > 90')
    df['low'] = df[['sstart', 'send']].min(axis=1)
    df['high'] =  df[['sstart', 'send']].max(axis=1)


    for toxin in toxin_genes:
        toxin_df = df.query('sseqid == @toxin')
        # print(toxin_df)

        # Now we only want to allow a contig (qseqid) to contribute it's best alignment
        toxin_df = df.loc[toxin_df.groupby('qseqid')['length'].idxmax()]
        toxin_df = toxin_df.sort_values(['low', 'high'], ascending=[True, False])
        # print(toxin_df)

        if toxin_df.shape[0] == 0:
            continue

        intervals = []
        identity_estimate = 0
        toxin_length = 0
        current_interval = (-1, -1)
        for low, high, pc_identity, slen in zip(toxin_df['low'], toxin_df['high'], toxin_df['pident'], toxin_df['slen']):
            if current_interval == (-1, -1):
                current_interval = (low, high)
                toxin_length = slen
                identity_estimate = pc_identity * ((high - low + 1) / toxin_length)
            elif low > current_interval[1]:
                # So this interval is to the right of previous
                intervals.append(current_interval)
                identity_estimate += pc_identity * ((high - low + 1) / toxin_length)
                current_interval = (low, high)
            elif high <= current_interval[1]:
                # Then this interval is completely within previous
                continue
            else:
                # This interval partially overlaps with previous
                identity_estimate += pc_identity * ((high - current_interval[1]) / toxin_length)
                current_interval = (current_interval[0], high)
            
        intervals.append(current_interval)
        
        covered_length = 0
        for start, end in intervals:
            covered_length += end + 1 - start
        identity_estimate = identity_estimate * (toxin_length / covered_length)

        if covered_length/toxin_length > 0.8 or identity_estimate > 97:
            toxin_dict[toxin] = {"presence": True, "percent_identity": identity_estimate, "length": covered_length, "gene_length": toxin_length}
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
    toxin_dict = find_toxin_genes_sum(catalogue, blast_output_tsv, toxin_dict)

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