#! /env/bin/python3
 
from argparse import ArgumentParser, SUPPRESS
from fpdf import FPDF
import logging
from pathlib import Path
import csv
import json

def generate_summary_report_args(parser):
    parser.add_argument('-s', '--samples_json', required=True,
                            help='Path to basic sample info TSV')
    parser.add_argument('-o', '--output_pdf', required=False, default="summary_report.pdf",
                            help='Path to output PDF')
    return parser

def generate_summary_report(samples_json_file: str, output_pdf:str):
    # Read in samples summary data in JSON
    samples_json = {}
    try:
        with open(samples_json_file) as file:
            samples_json = json.load(file)
    except IOError as e:
        logging.error(f"Error opening samples JSON {samples_json_file}")
        logging.error(e)
        exit(1)

    clusters = {}
    with open(samples_json["cluster_file"]) as cluster_file:
        cluster_reader = csv.DictReader(cluster_file, delimiter="\t")
        for line in cluster_reader:
            if line["cluster_number"] in clusters:
                clusters[line["cluster_number"]].append(line["id"])
            else:
                clusters[line["cluster_number"]] = [line["id"]]

    filtered_clusters = {"No Cluster" : []}
    for cluster_id, cluster in clusters.items():
        if len(cluster) < 3:
            filtered_clusters["No Cluster"].extend(cluster)
        else:
            filtered_clusters[cluster_id] = cluster
    clusters = filtered_clusters

    #Append logos
    WIDTH = 210
    #HEIGHT = 297
    pdf = FPDF()
    pdf.add_page()
    epw = pdf.w - 2*pdf.l_margin

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
    pdf.cell(w=0, h=5, txt="Sequencing QC summary", align = "L", ln=2)
    pdf.ln(1)

    # QC info
    pdf.set_text_color(0,0,0)
    pdf.cell(w=60, h=5, txt = "Sample Name", align = "C", border="TBL")
    pdf.cell(w=40, h=5, txt = "QC result", align = "C", border="TBR", ln=1)
    pdf.set_font("Helvetica", size=11)
    for sample, info in samples_json["samples"].items():
        qc_pass = "Pass"
        try:
            with open(info["qc_file"]) as file:
                qc_tsv_reader = csv.reader(file, delimiter="\t")
                for i, line in enumerate(qc_tsv_reader):
                    while("" in line):
                        line.remove("")
                    
                    if i > 0 and i <= 5:
                        if line[3] != "Pass":
                            qc_pass = "Fail"
                            break
            info["qc_pass"] = qc_pass
        except Exception as e:
            logging.error(f"Error opening qc_file '{info['qc_file']}' for sample {sample}")
            logging.error(e)
            exit(1)
        
        pdf.cell(w=60, h=5, txt = sample, align = "C", border="TBL")
        pdf.cell(w=40, h=5, txt = qc_pass, align = "C", border="TBR", ln=1)
    pdf.ln(5)

    # Cluster Info
    pdf.set_font("Helvetica", "B", size=12)
    pdf.set_text_color(0,0,200) #blue
    pdf.cell(w=0, h=5, txt="Sample clusters", align = "L", ln=2)
    pdf.ln(1)

    pdf.set_text_color(0,0,0)
    pdf.cell(w=epw/6, h=5, txt = "Cluster No", align = "C", border="TBL")
    pdf.cell(w=5*epw/6, h=5, txt = "Samples", align = "C", border="TBR", ln=1)
    pdf.set_font("Helvetica", size=11)
    line_width = 80

    for cluster_no, samples in clusters.items():
        cluster_str = [str(samples).strip("[]")]
        while len(cluster_str[-1]) > line_width:
            space = 0
            while cluster_str[-1][line_width - space] != ' ':
                space+=1
            last_str = cluster_str[-1][line_width - space:]
            cluster_str.append(last_str)
            cluster_str[-2] = cluster_str[-2][0:line_width - space]

        if len(cluster_str) > 1:
            for i, line in enumerate(cluster_str):
                if i == 0:
                    pdf.cell(w=epw/6, h=5, txt = cluster_no, align = "C", border="TL")
                    pdf.cell(w=5*epw/6, h=5, txt = line, align = "L", border="TR", ln=1)
                elif i < len(cluster_str) -1:
                    pdf.cell(w=epw/6, h=5, txt = "", align = "C", border="L")
                    pdf.cell(w=5*epw/6, h=5, txt = line, align = "L", border="R", ln=1)
                else:
                    pdf.cell(w=epw/6, h=5, txt = "", align = "C", border="BL")
                    pdf.cell(w=5*epw/6, h=5, txt = line, align = "L", border="BR", ln=1)
        else:
            pdf.cell(w=epw/6, h=5, txt = cluster_no, align = "C", border="TBL")
            pdf.cell(w=5*epw/6, h=5, txt = cluster_str[0], align = "L", border="TBR", ln=1)
    pdf.ln(1)
    pdf.set_font("Helvetica", "B", size=8)
    pdf.cell(7, 5, "N.B.")
    pdf.set_font("Helvetica", size=8)
    pdf.cell(WIDTH, 5, "Only clusters of 3 or more samples are considered. Clusters of 1 or 2 samples are listed under \"No Cluster\".")
    pdf.ln(5)

    pdf.set_font("Helvetica", "B", size=12)
    for cluster_no, samples in clusters.items():
        if cluster_no != "No Cluster":
            pdf.add_page()
            pdf.cell(w=0, h=5, txt=f"Cluster {cluster_no}", align = "L", ln=2)
            pdf.image(f"{samples_json['cluster_trees_dir']}/cluster_{cluster_no}_iqtree_scaled.png", w=3*WIDTH/4)
            pdf.ln(5)

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

