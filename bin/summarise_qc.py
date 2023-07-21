from argparse import ArgumentParser
import os
import json
import pandas as pd


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

def read_mixed_infection(dir):
    dfs = []
    for filename in os.listdir(dir):
        if filename.endswith('_mixed_infection_estimate.tsv'):
            file_path = os.path.join(dir, filename)

            df = pd.read_csv(file_path, sep='\t', index_col=False)
            dfs.append(df)
    
    df = pd.concat(dfs)
    return df

def read_mlst(dir):
    dfs = []
    for filename in os.listdir(dir):
        if filename.endswith('_name.tsv'):
            file_path = os.path.join(dir, filename)

            id = filename.replace('_name.tsv', '')
            df = pd.read_csv(file_path, sep='\t', index_col=False)
            df['id'] = id
            df = df[['id', 'MLST', 'Equivalent Ribotype', 'Collection Date', 'Source Hospital']]
            dfs.append(df)

    df = pd.concat(dfs)
    return df

def read_toxins(dir):
    rows = []
    for filename in os.listdir(dir):
        if filename.endswith('_toxin_coding_genes_report.json'):
            file_path = os.path.join(dir, filename)
            id = filename.replace('_toxin_coding_genes_report.json', '')

            with open(file_path, 'r') as json_file:
                data = json.load(json_file)
                # Create one long row
                row = {'id' : id}
                for gene, gene_data in data.items():
                    row[f"{gene}_presence"] = gene_data['presence']
                    row[f"{gene}_percent_identity"] = max(gene_data['percent_identity'], 0)
                    row[f"{gene}_length"] = max(gene_data['length'], 0)
                    row[f"{gene}_gene_length"] = gene_data['gene_length']
                rows.append(row)

    df = pd.DataFrame(rows)
    return df

def read_amr(dir):
    rows = []

    for filename in os.listdir(dir):
        if filename.endswith('_resistance_report.json'):
            file_path = os.path.join(dir, filename)
            id = filename.replace('_resistance_report.json', '')

            amr_df = pd.read_json(file_path).transpose()
            amr_df.index.names = ['Drug']
            amr_df = amr_df.reset_index()

            row = {'id' : id}
            for drug, res in zip(amr_df['Drug'], amr_df['resistance']):
                row[f"{drug}_resistance"] = res
            
            for drug, res_genes, missing_genes in zip(amr_df['Drug'], amr_df['evidence_resistance'], amr_df['evidence_sensitive']):
                for g in res_genes:
                    g = g.replace(' ', '')
                    row[f"{g} ({drug})"] = 1
                for g in missing_genes:
                    g = g.replace(' ', '')
                    row[f"{g} ({drug})"] = 0
            
            rows.append(row)
    
    return pd.DataFrame(rows).fillna(0)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-d', '--dir', required=True,
                        help='Directory containing qc tsv (should be reports folder)')
    parser.add_argument('-o', '--output', required=True,
                        help='Ouput file')
    args = parser.parse_args()
    dir = args.dir
    output_file = args.output

    qc_df = read_stadard_qc(dir)
    bracken_df = read_bracken(dir)
    mixed_infection_df = read_mixed_infection(dir)
    mlst_df = read_mlst(dir)
    toxin_df = read_toxins(dir)
    amr_df = read_amr(dir)

    df =  mlst_df.merge(qc_df, on='id', how='left') \
            .merge(mixed_infection_df, on='id', how='left') \
            .merge(bracken_df, on='id', how='left') \
            .merge(toxin_df, on='id', how='left') \
            .merge(amr_df, on='id', how='left')
    df.to_csv(output_file, index=False)