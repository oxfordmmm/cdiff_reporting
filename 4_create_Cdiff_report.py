# simple_table.py
 
from argparse import ArgumentParser, SUPPRESS
from fpdf import FPDF
import pandas as pd
import logging
from pathlib import Path

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
    #Sample name
    try:
        with open(sample_tsv) as h:
            sample = h.readlines()
        sample_name_str = ("").join(sample)
    except IOError as e:
        logging.error(f"Error opening sample TSV {sample_tsv}")
        logging.error(e)
        exit(1)

    #print(name)

    #Read in QC stats in TSV
    try:
        with open(qc_tsv) as g:
            lines = g.readlines()
        qc_tsv_str = ("").join(lines)
    except IOError as e:
        logging.error(f"Error opening QC TSV {qc_tsv}")
        logging.error(e)
        exit(1)
    #print(tsv)

    #Read in relatedness TSV
    try:
        with open(relatedness_tsv) as g:
            lines = g.readlines()
        relatedness_tsv_str = ("").join(lines)
    except IOError as e:
        logging.error(f"Error opening Relatedness TSV {relatedness_tsv}")
        logging.error(e)
        exit(1)

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
    pdf.set_font("Times", "B", size=15)
    pdf.ln(10)  # move 10 down
    pdf.cell(w=0, h=5, txt="Clostridioides difficile Sequence Analysis Report", align = "L", ln=2)
    #pdf.ln(5)

#Show sample details
    pdf.set_font("Times", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.ln(7)  # move 10 down
    pdf.cell(w=0, h=5, txt="Sample details", align = "L", ln=2)
    pdf.ln(5)

    pdf.set_font("Times", size=12)
    pdf.set_text_color(0,0,0)
    pdf.multi_cell(w=200, h=5, txt = sample_name_str)
    pdf.ln(5)

#Append QC stats
    pdf.set_font("Times", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Sequencing Quality Stats", align = "L", ln=2)
    
    pdf.set_font("Times", size=12)
    pdf.set_text_color(0,0,0)
    pdf.multi_cell(w=200, h=5, txt = qc_tsv_str)
    pdf.ln(5)

#Append AMR Profile
    pdf.set_font("Times", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Antimicrobial Resistance Profile", align = "L", ln=2)
    pdf.set_font("Times", size=12)
    pdf.ln(10)

    pdf.set_font("Times", size=12)
    pdf.set_text_color(0,0,0)
    #pdf.multi_cell(w=200, h=5, txt = str(data_df))
    epw = pdf.w - 2*pdf.l_margin
 
    row_height = pdf.font_size
    pdf.set_font("Times", "B", size=12)
    pdf.cell(epw/3, row_height, "Drug", border="TL", ln=0)
    pdf.cell(pdf.font_size * 3, row_height, "S/R", border="T", ln=0)
    pdf.cell((4*(epw/6)) - (pdf.font_size*3), row_height, "Evidence of Resistance", border="TR", ln=1)
    pdf.cell(epw, row_height, "Evidence of Sensitivity", border="BLR", ln=1)

    pdf.set_font("Times", size=12)
    for row in data_df.index:
        pdf.cell(epw/3, row_height, str(row), border="TL", ln=0)
        pdf.cell(pdf.font_size * 3, row_height, str(data_df['resistance'][row]), border="T", ln=0)
        pdf.cell((4*(epw/6)) - (pdf.font_size*3), row_height, str(data_df['evidence_resistance'][row]), border="TR", ln=1)
        pdf.cell(epw, row_height, str(data_df['evidence_sensitive'][row]), border="BLR", ln=1)
        
        
#Append Related samples section
    pdf.ln(5)
    pdf.set_font("Times", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Related samples", align = "L", ln=2)
    pdf.set_font("Times", size=14)
    pdf.ln(5)

    pdf.set_font("Times", size=12)
    pdf.set_text_color(0,0,0)
    pdf.multi_cell(w=200, h=5, txt = relatedness_tsv_str)
    pdf.ln(5)

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

