#! /env/bin/python3
 
from argparse import ArgumentParser, SUPPRESS
from fpdf import FPDF
import pandas as pd
import logging
from pathlib import Path
import csv

def generate_summary_report_args(parser):
    parser.add_argument('-s', '--samples_json', required=True,
                            help='Path to basic sample info TSV')
    parser.add_argument('-o', '--output_pdf', required=False, default="summary_report.pdf",
                            help='Path to output PDF')
    return parser

def generate_summary_report(samples_json: str, output_pdf:str):
    # Read in samples summary data in JSON
    try:
        pd.set_option('display.max_colwidth', None)
        amr_df = pd.read_json(samples_json)
    except IOError as e:
        logging.error(f"Error opening samples JSON {samples_json}")
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
    pdf.cell(w=0, h=5, txt="Clostridioides difficile Cluster Summary Report", align = "L", ln=2)
    pdf.ln(5)

    #Show sample details
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Samples summary", align = "L", ln=2)
    pdf.ln(1)

    #Save to PDF
    try:
        pdf.output(output_pdf)
    except Exception as e:
        logging.error(f"Error writing report to {output_pdf}")
        logging.error(e)
        exit(1)

if __name__ == '__main__':
    parser = ArgumentParser(description='Produce a PDF showing a summary of Cdiff sample clusters.')
    parser = generate_summary_report_args(parser)
    args = parser.parse_args()

    logging.basicConfig(handlers=[
            logging.StreamHandler()],
        format='%(asctime)s - %(levelname)s - %(message)s', 
        level=logging.INFO)

    generate_summary_report(args.samples_json, args.output_pdf)

