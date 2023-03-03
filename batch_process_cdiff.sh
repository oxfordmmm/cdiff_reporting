SOURCE_DIR="${PWD}/sample_data"
OUTPUT_DIR="${PWD}/sample_data/report"
CGMLST_DIR="${PWD}/sample_data/cgmlst"
DATA_DIR="${PWD}/data"

if [ -f "$OUTPUT_DIR" ]
then 
    echo "$OUTPUT_DIR is not empty"
    exit 1
fi

for path in $(ls ${SOURCE_DIR}/assemblies/*_R_contigs.fa) 
do 
    basename $path | perl -p -e 's/\_contigs\.fa//g' >> "${OUTPUT_DIR}/processed_genomes.tsv"
done

for genome in $(cat ${OUTPUT_DIR}/processed_genomes.tsv)
do
    #rm QC_summary_table.tsv
    #sample details eg name, mlst
    perl bin/get_sample_name_MLST.pl "${SOURCE_DIR}/mlst/${genome}_ST.tsv" "${OUTPUT_DIR}/${genome}_name.tsv" "${DATA_DIR}/LookUpTable_ST_RT.tsv"
    #qc
    perl bin/parseQCoutputs.pl "${SOURCE_DIR}/raw_fastqc_single/${genome}_raw_reads_fastqc/${genome}.txt" "${SOURCE_DIR}/quast/${genome}_Quastreport.tsv" "${OUTPUT_DIR}/${genome}_QC_summary_table.tsv"
    #amr profile
    python3 bin/AMR_process.py -c data/AMR_catalogue.json \
    -s data/AMR_catalogue_schema.json \
    -b "${SOURCE_DIR}/cdiff_blastn/cdiffamr-${genome}_blastn.tsv" \
    -f "${SOURCE_DIR}/amr_cdiff_pointmuts/${genome}_forpointmuts.tsv"  \
    -o "${OUTPUT_DIR}/${genome}_resistance_report.json"
    #relatedness
    python3 bin/compareProfilesSingleExclude.py -f "${CGMLST_DIR}" \
    -i "${SOURCE_DIR}/cgmlst/${genome}_cgmlst.json" \
    -d 20 \
    -n 3 \
    -o "${OUTPUT_DIR}/${genome}_cgmlst_comparisons.tsv"
    #final report
    python3 bin/AMR_report.py -s "${OUTPUT_DIR}/${genome}_name.tsv" \
    -q "${OUTPUT_DIR}/${genome}_QC_summary_table.tsv" -a "${OUTPUT_DIR}/${genome}_resistance_report.json" \
    -r "${OUTPUT_DIR}/${genome}_cgmlst_comparisons.tsv" \
    -o "${OUTPUT_DIR}/${genome}_out_report.pdf"
done