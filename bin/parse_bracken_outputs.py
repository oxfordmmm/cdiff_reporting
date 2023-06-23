#! /env/bin/python3
from argparse import ArgumentParser
import os
import sys
import pandas as pd

phage_list = ['Clostridium phage phiCD6356', 'Colneyvirus CD27', 'Colneyvirus CDKM15',
              'Colneyvirus CDKM9', 'Colneyvirus MMP02', 'Leicestervirus CD111', 'Leicestervirus CD146',
              'Leicestervirus CD382', 'Lubbockvirus CD119', 'Lubbockvirus CDHM19', 'Sherbrookevirus CD4811',
              'Sherbrookevirus CD506', 'Sherbrookevirus CDHM13', 'Sherbrookevirus CDHM14', 'Sherbrookevirus MMP04',
              'Yongloolinvirus CDMH1', 'Yongloolinvirus MMP01']

plausible_misclassification = ['Enterocloster bolteae', '[Clostridium] scindens',
        '[Clostridium] innocuum', '[Ruminococcus] gnavus',
        'Roseburia intestinalis', 'Blautia obeum', 'Enterocloster clostridioformis',
        'Clostridium perfringens',
        'Clostridioides sp. ES-S-0107-01',
        'Eubacterium maltosivorans',
        'Anaerostipes caccae',
        'Romboutsia ilealis',
        'Clostridioides sp. ES-S-0173-01',
        'Peptacetobacter hiranonis',
        'Flavobacterium sp. N502536',
        'Dorea longicatena',
        'Flavonifractor plautii',
        'Clostridium cadaveris',
        'Qiania dongpingensis',
        'Blautia sp. NBRC 113351',
        'Blautia producta',
        'Blautia obeum',
        'Blautia argi',
        'Coprococcus comes',
        'Paeniclostridium sordellii',
        'Faecalibacterium sp. IP-3-29',
        'Thomasclavelia ramosa',
        'Blautia wexlerae'
]   

def flatten(l):
    return [item for sublist in l for item in sublist] 

def get_bracken_top_hits(bracken_df, num):
    df = bracken_df.iloc[:num]
    return [(name, round(100*pc,2)) for name, pc in zip(df['name'], df['fraction_total_reads'])]

def get_group_sums(bracken_df):
    cdiff_pc = 100 * bracken_df.loc[bracken_df['name'] == 'Clostridioides difficile', 'fraction_total_reads'].sum()

    phage_mask = bracken_df['name'].isin(phage_list)
    phage_pc = 100 * bracken_df.loc[phage_mask, 'fraction_total_reads'].sum()

    misclassification_mask = bracken_df['name'].isin(plausible_misclassification)
    misclassification_pc = 100 * bracken_df.loc[misclassification_mask, 'fraction_total_reads'].sum()

    other_pc = 100 - cdiff_pc - phage_pc - misclassification_pc
    return [cdiff_pc, phage_pc, misclassification_pc, other_pc]

def write_qc(group_sums, top_hits, outfile, num):
    with open(outfile, 'w') as file:
        group_cols = ['cdiff_pc', 'phage_pc', 'misclassification_pc', 'other_pc']
        top_hit_cols = flatten([[f"bracken_{str(i)}", f"pc_{str(i)}"]  for i in range(1,num+1)])

        group_values = [str(round(pc, 2)) for pc in group_sums]
        top_hit_values = flatten([ [species, str(pc)] for species, pc in top_hits])

        file.write('\t'.join(group_cols + top_hit_cols) + '\n')
        file.write('\t'.join(group_values + top_hit_values) + '\n')


if __name__ == '__main__':
    parser = ArgumentParser(description='Summarise bracken outputs to tsv.')
    parser.add_argument('-b', '--bracken_output', required=True)
    parser.add_argument('-o', '--outfile', required=True)
    parser.add_argument('-n', '--number', required=False, default=3)
    args = parser.parse_args()
    bracken_file = args.bracken_output
    num = int(args.number)

    if not os.path.exists(bracken_file):
        print("Error: bracken file doesn't exist.")
        grouped_pcs = [0,0,0,0]
        top_hits = [("", 0) for i in range(num)]
        write_qc(grouped_pcs, top_hits, args.outfile, num)
        sys.exit()

    try:
        bracken_df = pd.read_csv(bracken_file, sep='\t').sort_values('fraction_total_reads', ascending=False, ignore_index=True)
        top_hits = get_bracken_top_hits(bracken_df, num)
        grouped_pcs = get_group_sums(bracken_df)

        write_qc(grouped_pcs, top_hits, args.outfile, num)
    except Exception as e:
        print("Error occurred while processing bracken output:", str(e))


