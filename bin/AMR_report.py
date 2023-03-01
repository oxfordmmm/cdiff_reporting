# simple_table.py
 
from argparse import ArgumentParser, SUPPRESS
from fpdf import FPDF
import pandas as pd
import logging
from pathlib import Path
import csv

def amr_report_args(parser):
    parser.add_argument('-s', '--sample_tsv', required=True,
                            help='Path to basic sample info TSV')
    parser.add_argument('-q', '--qc_tsv', required=True,
                            help='Path to QC summary TSV')
    parser.add_argument('-a', '--amr_json', required=True,
                            help='Path to AMR JSON')
    parser.add_argument('-r', '--relatedness_tsv', required=True,
                            help='Path to relatedness TSV')
    parser.add_argument('-o', '--output_pdf', required=False, default="sample_report.pdf",
                            help='Path to output PDF')
    return parser

def amr_report(sample_tsv: str, qc_tsv: str, amr_json:str, relatedness_tsv:str, output_pdf:str):
    #Read in AMR data in JSON
    try:
        pd.set_option('display.max_colwidth', None)
        data_df = pd.read_json(amr_json).transpose()
        data_df.index.names = ['Drug']
    except IOError as e:
        logging.error(f"Error opening AMR JSON {amr_json}")
        logging.error(e)
        exit(1)
    #print(data_df)

    #Append logos
    WIDTH = 210
    #HEIGHT = 297
    pdf = FPDF()
    pdf.add_page()
    pdf.image('data/oxford.png', x=0, y=0, w=WIDTH/2-10)
    pdf.image('data/leeds.jpg', x=120, y=WIDTH/9.5, w=80)
    pdf.ln(30)  # move 30 down

    #Show poject title
    pdf.set_font("Helvetica", "B", size=17)
    pdf.ln(10)  # move 10 down
    pdf.cell(w=0, h=5, txt="Clostridioides difficile Sequence Analysis Report", align = "L", ln=2)
    pdf.ln(5)

    #Show sample details
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Sample details", align = "L", ln=2)
    pdf.ln(1)

    #pdf.set_font("Helvetica", size=12)
    pdf.set_text_color(0,0,0)

    try:
        with open(sample_tsv) as file:
            sample_name_tsv = csv.reader(file, delimiter="\t")
            
            for i, line in enumerate(sample_name_tsv):
                while("" in line):
                    line.remove("")
                
                if i == 0:
                    pdf.set_font("Helvetica", "B", size=12)
                    pdf.cell(w=40, h=5, txt = line[0], border="TBL")
                    pdf.cell(w=60, h=5, txt = line[1], border="TB")
                    pdf.cell(w=60, h=5, txt = line[2], border="TB", ln=1)
            
                else:
                    pdf.set_font("Helvetica", size=12)
                    pdf.cell(w=40, h=5, txt = line[0], border="TBL")
                    pdf.cell(w=60, h=5, txt = line[1], border="TB")
                    pdf.cell(w=60, h=5, txt = line[2], border="TB", ln=1)
                    
            pdf.ln(10)
    except Exception as e:
        logging.error(f"Error opening sample TSV {sample_tsv}")
        logging.error(e)
        exit(1)

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
                
                pdf.cell(w=60, h=5, txt = line[0], border="TBL")
                pdf.cell(w=40, h=5, txt = line[1], border="TB")
                pdf.cell(w=30, h=5, txt = line[2], border="TB")
                pdf.cell(w=30, h=5, txt = line[3], border="TBR", ln=1)
                pdf.set_font("Helvetica", size=12)
    except IOError as e:
        logging.error(f"Error opening QC TSV {qc_tsv}")
        logging.error(e)
        exit(1)

    pdf.ln(10)

    #Append AMR Profile
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Antimicrobial Resistance Profile", align = "L", ln=1)
    pdf.set_font("Helvetica", size=12)

    pdf.set_font("Helvetica", size=12)
    pdf.set_text_color(0,0,0)
    epw = pdf.w - 2*pdf.l_margin
 
    row_height = pdf.font_size
    pdf.set_font("Helvetica", "B", size=12)
    pdf.cell(epw/3, row_height, "Drug", border="TL", ln=0)
    pdf.cell(pdf.font_size * 3, row_height, "S/R", border="T", ln=0)
    pdf.cell((4*(epw/6)) - (pdf.font_size*3), row_height, "Evidence of Resistance", border="TR", ln=1)
    pdf.cell(epw, row_height, "Catalogue features not found", border="BLR", ln=1)

    pdf.set_font("Helvetica", size=12)
    for row in data_df.index:
        pdf.cell(epw/3, row_height, str(row), border="TL", ln=0)
        pdf.cell(pdf.font_size * 3, row_height, str(data_df['resistance'][row]), border="T", ln=0)
        pdf.cell((4*(epw/6)) - (pdf.font_size*3), row_height, str(data_df['evidence_resistance'][row]), border="TR", ln=1)
        pdf.cell(epw, row_height, str(data_df['evidence_sensitive'][row]), border="BLR", ln=1)

    #Append Related samples section
    pdf.ln(10)
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Related samples", align = "L", ln=2)

    #Read in relatedness TSV
    pdf.set_text_color(0,0,0) #black
    try:
        with open(relatedness_tsv) as file:
            relatedness_tsv_reader = csv.reader(file, delimiter="\t")
            for line in relatedness_tsv_reader:
                print(line)
                while("" in line):
                    line.remove("")

                if len(line) < 4:
                    print(f"sample {line[0]} outside cutoff")
                    continue
                pdf.cell(w=90 , h=5, txt = line[0], border="TBL")
                pdf.cell(w=40 , h=5, txt = line[1], border="TB")
                pdf.cell(w=30 , h=5, txt = line[2], border="TB")
                pdf.cell(w=30 , h=5, txt = line[3], border="TBR", ln=1)
                pdf.set_font("Helvetica", size=12)
    except IOError as e:
        logging.error(f"Error opening Relatedness TSV {relatedness_tsv}")
        logging.error(e)
        exit(1)

    #Save to PDF
    try:
        pdf.output(output_pdf)
    except Exception as e:
        logging.error(f"Error writing report to {output_pdf}")
        logging.error(e)
        exit(1)

if __name__ == '__main__':
    parser = ArgumentParser(description='Produce a PDF showing basic, QC, AMR and relatedness info for Cdiff samples.')
    parser = amr_report_args(parser)
    args = parser.parse_args()

    logging.basicConfig(handlers=[
            logging.StreamHandler()],
        format='%(asctime)s - %(levelname)s - %(message)s', 
        level=logging.INFO)

    amr_report(args.sample_tsv, args.qc_tsv, args.amr_json, args.relatedness_tsv, args.output_pdf)

