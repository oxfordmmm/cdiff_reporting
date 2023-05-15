# simple_table.py
 
from argparse import ArgumentParser, SUPPRESS
from fpdf import FPDF
import pandas as pd
import logging
from pathlib import Path
import csv

def generate_individual_report_args(parser):
    parser.add_argument('-s', '--sample_tsv', required=True,
                            help='Path to basic sample info TSV')
    parser.add_argument('-q', '--qc_tsv', required=True,
                            help='Path to QC summary TSV')
    parser.add_argument('-a', '--amr_json', required=True,
                            help='Path to AMR JSON')
    parser.add_argument('-t', '--toxin_json', required=True,
                            help='Path to toxin gene presence JSON')
    parser.add_argument('-r', '--relatedness_tsv', required=True,
                            help='Path to relatedness TSV')
    parser.add_argument('-o', '--output_pdf', required=False, default="sample_report.pdf",
                            help='Path to output PDF')
    return parser

def generate_individual_report(sample_tsv: str, qc_tsv: str, amr_json:str, toxin_json:str, relatedness_tsv:str, output_pdf:str):
    # Read in AMR data in JSON
    try:
        pd.set_option('display.max_colwidth', None)
        amr_df = pd.read_json(amr_json).transpose()
        amr_df.index.names = ['Drug']
    except IOError as e:
        logging.error(f"Error opening AMR JSON {amr_json}")
        logging.error(e)
        exit(1)
    #print(amr_df)

    # Read in toxin coding genes in JSON
    try:
        toxin_df = pd.read_json(toxin_json, typ="series").transpose()
    except IOError as e:
        logging.error(f"Error opening AMR JSON {toxin_json}")
        logging.error(e)
        exit(1)

    #Append logos
    WIDTH = 210
    #HEIGHT = 297
    pdf = FPDF()
    pdf.add_page()
    pdf.image('data/oxford.png', x=0, y=0, w=WIDTH/2-10)
    pdf.image('data/leeds.jpg', x=120, y=WIDTH/9.5, w=80)
    pdf.ln(30)  # move 30 down

    #Show poject title
    pdf.set_font("Helvetica", "B", size=16)
    pdf.ln(10)  # move 10 down
    pdf.cell(w=0, h=5, txt="Clostridioides difficile Sequence Analysis Report", align = "L", ln=2)
    pdf.ln(5)

    #Show sample details
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Sample details", align = "L", ln=2)
    pdf.ln(1)

    #pdf.set_font("Helvetica", size=11)
    pdf.set_text_color(0,0,0)

    try:
        with open(sample_tsv) as file:
            #sample_name_tsv = csv.reader(file, delimiter="\t")
            sample_name_tsv_df = pd.read_csv(file, sep = '\t')

            # Specimen Identifier
            pdf.set_font("Helvetica", "B", size=11)
            pdf.cell(w=40, h=5, txt = sample_name_tsv_df.columns[1], border="TBL")
            pdf.set_font("Helvetica", size=10)
            pdf.cell(w=60, h=5, txt = str(sample_name_tsv_df[sample_name_tsv_df.columns[1]][0]), border="TBR")
            # Source Hospital
            pdf.set_font("Helvetica", "B", size=11)
            pdf.cell(w=40, h=5, txt = sample_name_tsv_df.columns[5], border="TB")
            pdf.set_font("Helvetica", size=10)
            pdf.cell(w=50, h=5, txt = str(sample_name_tsv_df[sample_name_tsv_df.columns[5]][0]), border="TBR", ln=1)
    
            # Report date
            pdf.set_font("Helvetica", "B", size=11)
            pdf.cell(w=40, h=5, txt = sample_name_tsv_df.columns[0], border="TBL")
            pdf.set_font("Helvetica", size=10)
            pdf.cell(w=60, h=5, txt = str(sample_name_tsv_df[sample_name_tsv_df.columns[0]][0]), border="TBR")
            # Collection Date
            pdf.set_font("Helvetica", "B", size=11)
            pdf.cell(w=40, h=5, txt = sample_name_tsv_df.columns[4], border="TB")
            pdf.set_font("Helvetica", size=10)
            pdf.cell(w=50, h=5, txt = str(sample_name_tsv_df[sample_name_tsv_df.columns[4]][0]), border="TBR", ln=1)

            # MLST
            pdf.set_font("Helvetica", "B", size=11)
            pdf.cell(w=15, h=5, txt = sample_name_tsv_df.columns[2], border="TBL")
            pdf.set_font("Helvetica", size=10)
            pdf.cell(w=7, h=5, txt = str(sample_name_tsv_df[sample_name_tsv_df.columns[2]][0]), border="TBR")
            # Ribotype
            pdf.set_font("Helvetica", "B", size=11)
            pdf.cell(w=50, h=5, txt = sample_name_tsv_df.columns[3], border="TB")
            pdf.set_font("Helvetica", size=9)
            pdf.cell(w=28, h=5, txt = str(sample_name_tsv_df[sample_name_tsv_df.columns[3]][0]), border="TBR")
            # Pipeline version
            pdf.set_font("Helvetica", "B", size=11)
            pdf.cell(w=40, h=5, txt = "Pipeline Version", border="TB")
            pdf.set_font("Helvetica", size=10)
            pdf.cell(w=50, h=5, txt = "v0.0.1", border="TBR", ln=1)
    except Exception as e:
        logging.error(f"Error opening sample TSV {sample_tsv}")
        logging.error(e)
        exit(1)
    pdf.ln(5)

    #Append QC stats
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Sequencing Quality Stats", align = "L", ln=2)
    
    pdf.set_text_color(0,0,0)
    #Read in QC stats in TSV
    try:
        with open(qc_tsv) as file:
            qc_tsv_reader = csv.reader(file, delimiter="\t")
            for i, line in enumerate(qc_tsv_reader):
                while("" in line):
                    line.remove("")
                
                pdf.cell(w=60, h=5, txt = line[0], align = "L", border="TBL")
                pdf.cell(w=40, h=5, txt = line[1], align = "C", border="TB")
                pdf.cell(w=30, h=5, txt = line[2], align = "C", border="TB")
                pdf.cell(w=30, h=5, txt = line[3], align = "C", border="TBR", ln=1)
                pdf.set_font("Helvetica", size=11)
    except IOError as e:
        logging.error(f"Error opening QC TSV {qc_tsv}")
        logging.error(e)
        exit(1)

    pdf.ln(5)

    # Toxin coding genes P/A
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Toxin-coding Genes", align = "L", ln=1)
    pdf.set_font("Helvetica", size=11)

    pdf.set_font("Helvetica", size=11)
    pdf.set_text_color(0,0,0)
    epw = pdf.w - 2*pdf.l_margin

    pdf.set_font("Helvetica", "B", size=11)
    pdf.cell(30, 5, "Toxin Gene", border="TBL", ln=0)
    pdf.cell(20, 5, "Presence", border="TB", align = "C", ln=0)
    pdf.cell(epw/6, 5, "Identity %", border="TB", align = "C", ln=0)
    pdf.cell(epw/6, 5, "Match Length", border="TB", align = "C", ln=0)
    pdf.cell(epw/6, 5, "Match %", border="TB", align = "C", ln=0)
    pdf.cell(epw/6, 5, "Gene Length", border="TBR", align = "C", ln=1)

    pdf.set_font("Helvetica", size=11)
    for row in toxin_df.index:
        pdf.cell(30, 5, row, border="TBL", ln=0)
        if toxin_df[row]["presence"]:
            pdf.cell(20, 5, "Present", border="TB", align = "C", ln=0)
            pdf.cell(epw/6, 5, str(toxin_df[row]["percent_identity"]), border="TB", align = "C", ln=0)
            pdf.cell(epw/6, 5, str(toxin_df[row]["length"]), border="TB", align = "C", ln=0)
            pdf.cell(epw/6, 5, str(toxin_df[row]["length"] / toxin_df[row]["gene_length"]), border="TBR", align = "C", ln=1)
            pdf.cell(epw/6, 5, str(toxin_df[row]["gene_length"]), border="TBR", align = "C", ln=1)
        else:
            pdf.cell(20, 5, "Not Found", border="TB", align = "C", ln=0)
            pdf.cell(epw/6, 5, "N/A", border="TB", align = "C", ln=0)
            pdf.cell(epw/6, 5, "N/A", border="TB", align = "C", ln=0)
            pdf.cell(epw/6, 5, "N/A", border="TB", align = "C", ln=0)
            pdf.cell(epw/6, 5, str(toxin_df[row]["gene_length"]), border="TBR", align = "C", ln=1)

    drug_str = ["Catalogue Features not found: ", ""]
    for row in amr_df.index:
        if len(amr_df['evidence_sensitive'][row]) > 0:
            drug_str[1] += f"{row} = {amr_df['evidence_sensitive'][row]}; ".translate( {ord(i): None for i in "[]"} )
        else:
            drug_str[1] += f"{row} = None; "

    pdf.ln(5)

    #Append AMR Profile
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Antimicrobial Resistance Profile", align = "L", ln=1)
    pdf.set_font("Helvetica", size=11)

    pdf.set_font("Helvetica", size=11)
    pdf.set_text_color(0,0,0)

    #row_height = pdf.font_size
    pdf.set_font("Helvetica", "B", size=11)
    pdf.cell(epw/3, 5, "Drug", border="TBL", ln=0)
    pdf.cell(pdf.font_size * 3, 5, "S/R", border="TB", align = "C", ln=0)
    pdf.cell((4*(epw/6)) - (pdf.font_size*3), 5, "Evidence of Resistance", border="TBR", align = "C", ln=1)
    #pdf.cell(epw, row_height, "Catalogue features not found", border="BLR", ln=1)

    pdf.set_font("Helvetica", size=11)
    for row in amr_df.index:
        pdf.cell(epw/3, 5, str(row), border="TBL", ln=0)
        pdf.cell(pdf.font_size * 3, 5, str(amr_df['resistance'][row]), border="TB", align = "C", ln=0)
        pdf.cell((4*(epw/6)) - (pdf.font_size*3), 5, str(amr_df['evidence_resistance'][row]).strip("[]"), border="TBR", align = "C", ln=1)
        #pdf.cell(epw, row_height, str(amr_df['evidence_sensitive'][row]), border="BLR", ln=1)

    drug_str = ["Catalogue Features not found: ", ""]
    for row in amr_df.index:
        if len(amr_df['evidence_sensitive'][row]) > 0:
            drug_str[1] += f"{row} = {amr_df['evidence_sensitive'][row]}; ".translate( {ord(i): None for i in "[]"} )
        else:
            drug_str[1] += f"{row} = None; "
    
    # Footnote features not found for sensitivity recording
    first_line = 1
    line_width = 140
    while len(drug_str[-1]) > line_width - (first_line * len(drug_str[0])):
        space = 0
        while drug_str[-1][line_width-(first_line * len(drug_str[0])) - space] != ' ':
            space+=1
        last_str = drug_str[-1][line_width-(first_line * len(drug_str[0])) - space:]
        drug_str.append(last_str)
        drug_str[-2] = drug_str[-2][0:line_width-(first_line * len(drug_str[0])) - space]
        first_line = 0
    
    pdf.set_font("Helvetica", "B", size=8)
    pdf.cell(45, 5, drug_str.pop(0))
    pdf.set_font("Helvetica", size=8)
    while len(drug_str) > 0:
        pdf.cell(WIDTH, 5, drug_str.pop(0), ln=1)
    pdf.ln(5)

    #Append Related samples section
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Related samples", align = "L", ln=2)

    #Read in relatedness TSV
    pdf.set_text_color(0,0,0) #black
    try:
        with open(relatedness_tsv) as file:
            relatedness_tsv_reader = csv.reader(file, delimiter="\t")
            for line in relatedness_tsv_reader:
                while("" in line):
                    line.remove("")

                if len(line) < 4:
                    print(f"sample {line[0]} outside cutoff")
                    continue
                pdf.cell(w=90 , h=5, txt = line[0], border="TBL")
                pdf.cell(w=40 , h=5, txt = line[1], align="C", border="TB")
                pdf.cell(w=30 , h=5, txt = line[2], align="C", border="TB")
                pdf.cell(w=30 , h=5, txt = line[3], align="C", border="TBR", ln=1)
                pdf.set_font("Helvetica", size=11)
    except IOError as e:
        logging.error(f"Error opening Relatedness TSV {relatedness_tsv}")
        logging.error(e)
        exit(1)
    pdf.ln(5)

    #Append Related samples section
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    #pdf.cell(w=0, h=5, txt="Cluster Report", align = "L", ln=1)

    #Save to PDF
    try:
        pdf.output(output_pdf)
    except Exception as e:
        logging.error(f"Error writing report to {output_pdf}")
        logging.error(e)
        exit(1)

if __name__ == '__main__':
    parser = ArgumentParser(description='Produce a PDF showing basic, QC, AMR and relatedness info for Cdiff samples.')
    parser = generate_individual_report_args(parser)
    args = parser.parse_args()

    logging.basicConfig(handlers=[
            logging.StreamHandler()],
        format='%(asctime)s - %(levelname)s - %(message)s', 
        level=logging.INFO)

    generate_individual_report(args.sample_tsv, args.qc_tsv, args.amr_json, args.toxin_json, args.relatedness_tsv, args.output_pdf)

