# simple_table.py
 
from fpdf import FPDF
import json
import pandas as pd
import numpy as np


def simple_table(spacing=1):

#Sample name
    with open('sample_name.tsv') as h:
        sample = h.readlines()
    name = ("").join(sample)
    print(name)

#Read in QC stats in TSV
    with open('QC_summary_table.tsv') as g:
        lines = g.readlines()
    tsv = ("").join(lines)
    print(tsv)
   
#Read in QC stats in JSON
    #data = pd.read_json('ERR3274281_resistance.json')
    pd.set_option('display.max_colwidth', None)
    data_df = pd.read_json('ERR3274281_resistance.json').transpose()
    data_df.index.names = ['Drug']
    print(data_df)

#Append logos
    WIDTH = 210
    #HEIGHT = 297
    pdf = FPDF()
    pdf.add_page()
    pdf.image('oxford.png', x=0, y=0, w=WIDTH/2-10)
    pdf.image('leeds.jpg', x=120, y=WIDTH/9.5, w=80)
    pdf.ln(30)  # move 30 down

#Show poject title
    pdf.set_font("Times", "B", size=17)
    pdf.ln(10)  # move 10 down
    pdf.cell(w=0, h=5, txt="Clostridioides difficile Sequence Analysis Report", align = "L", ln=2)
    pdf.ln(5)

#Show sample details
    pdf.set_font("Times", "B", size=14)
    pdf.set_text_color(0,0,200) #blue
    pdf.ln(15)  # move 10 down
    pdf.cell(w=0, h=5, txt="Sample details", align = "L", ln=2)
    pdf.ln(10)

    pdf.set_font("Times", size=14)
    pdf.set_text_color(0,0,0)
    pdf.multi_cell(w=200, h=5, txt = name)
    pdf.ln(10)

#Append QC stats
    pdf.set_font("Times", "B", size=14)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Sequencing Quality Stats", align = "L", ln=2)
    
    pdf.set_font("Times", size=14)
    pdf.set_text_color(0,0,0)
    pdf.multi_cell(w=200, h=5, txt = tsv)
    pdf.ln(10)

#Append AMR Profile
    pdf.set_font("Times", "B", size=14)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Antimicrobial Resistance Profile", align = "L", ln=2)
    pdf.set_font("Times", size=14)
    pdf.ln(10)

    pdf.set_font("Times", size=14)
    pdf.set_text_color(0,0,0)
    #pdf.multi_cell(w=200, h=5, txt = str(data_df))
    epw = pdf.w - 2*pdf.l_margin
 
    row_height = pdf.font_size
    pdf.set_font("Times", "B", size=14)
    pdf.cell(epw/3, row_height, "Drug", border='TL', ln=0)
    pdf.cell(pdf.font_size * 3, row_height, "S/R", border='T', ln=0)
    pdf.cell((4*(epw/6)) - (pdf.font_size*3), row_height, "Evidence of Resistance", border='TR', ln=1)
    pdf.cell(epw, row_height, "Evidence of Sensitivty", border='BLR', ln=1)

    pdf.set_font("Times", size=14)
    for row in data_df.index:
        pdf.cell(epw/3, row_height, str(row), border='TL', ln=0)
        pdf.cell(pdf.font_size * 3, row_height, str(data_df['resistance'][row]), border='T', ln=0)
        pdf.cell((4*(epw/6)) - (pdf.font_size*3), row_height, str(data_df['evidence_resistance'][row]), border='TR', ln=1)
        pdf.cell(epw, row_height, str(data_df['evidence_sensitive'][row]), border='BLR', ln=1)
    
#Save to PDF
    
    pdf.output('simple_table.pdf')
 
if __name__ == '__main__':
    simple_table()

