#! /bin/bash

unset $SOURCE_DIR
unset $DATA_DIR
unset $CGMLST_DIR
unset $OUTPUT_DIR
unset $MAPPING_DIR

while getopts "s:d:c:o:m:" flag
do
    case "${flag}" in
        s)
            if ! [[ -d "$OPTARG" ]]; then
                echo "Please enter a valid source directory."
                exit 1
            fi
            SOURCE_DIR=${OPTARG}
            ;;
        m)
            if ! [[ -d "$OPTARG" ]]; then
                echo "Please enter a valid mapping directory."
                exit 1
            fi
            MAPPING_DIR=${OPTARG}
            ;;
        d) 
            if ! [[ -d "$OPTARG" ]]; then
                echo "Please enter a valid data directory."
                exit 1
            fi
            DATA_DIR=${OPTARG}
            ;;
        c) 
            if ! [[ -d "$OPTARG" ]]; then
                echo "Please enter a valid cgmlst directory."
                exit 1
            fi
            CGMLST_DIR=${OPTARG}
            ;;
        o) 
            if ! [[ -d "$OPTARG" ]]; then
                mkdir -p ${OPTARG}/
            fi

            if [[ "$(ls -A $OPTARG)" ]]; then
                echo "The output directory is not empty."
                exit 1
            fi
            OUTPUT_DIR=${OPTARG}
            ;;
        
        :) usage 1 "-$OPTARG requires an argument" ;;
        ?) usage 1 "Unknown option '$opt'" ;;
    esac
done

if [[ -z ${SOURCE_DIR} ]]; then
    echo "Please enter a source directory"
    exit 1
elif [[ -z ${DATA_DIR} ]]; then
    echo "Please enter a data directory"
    exit 1
elif [[ -z ${CGMLST_DIR} ]]; then
    echo "Please enter a cgmlst directory"
    exit 1
elif [[ -z ${OUTPUT_DIR} ]]; then
    echo "Please enter a output directory"
    exit 1
fi

#SOURCE_DIR="${PWD}/sample_data"
#DATA_DIR="${PWD}/data"
#CGMLST_DIR="${PWD}/sample_data/cgmlst"
#OUTPUT_DIR="${PWD}/sample_data/report"

#mkdir ${OUTPUT_DIR}/

for path in $(ls ${SOURCE_DIR}/assemblies/*_contigs.fa) 
do 
    basename $path | perl -p -e 's/\_contigs\.fa//g' >> "${OUTPUT_DIR}/processed_genomes.tsv"
done

for genome in $(cat ${OUTPUT_DIR}/processed_genomes.tsv)
do
    #rm QC_summary_table.tsv
    #sample details eg name, mlst
    perl bin/get_sample_name_MLST.pl "${SOURCE_DIR}/mlst/${genome}_ST.tsv" "${OUTPUT_DIR}/${genome}_name.tsv" "${DATA_DIR}/LookUpTable_ST_RT.tsv"

    #qc
    python3 bin/parseQCoutputs.py -f "${SOURCE_DIR}/raw_fastqc_single/${genome}_raw_reads_fastqc/${genome}.txt" -q "${SOURCE_DIR}/quast/${genome}_Quastreport.tsv" -o "${OUTPUT_DIR}/${genome}_QC_summary_table.tsv"
    python3 bin/parse_bracken_outputs.py -b "${SOURCE_DIR}/kraken2/${genome}_bracken_report.tsv" -o "${OUTPUT_DIR}/${genome}_bracken_summary_table.tsv"

    #amr profile
    python3 bin/AMR_process.py -c data/AMR_catalogue.json \
    -s data/AMR_catalogue_schema.json \
    -b "${SOURCE_DIR}/cdiff_blastn/cdiffamr-${genome}_blastn.tsv" \
    -f "${SOURCE_DIR}/amr_cdiff_pointmuts/${genome}_forpointmuts.tsv"  \
    -o "${OUTPUT_DIR}/${genome}_resistance_report.json"

    #toxin coding genes
    python3 bin/toxin_coding_genes_process.py -c data/toxin_coding_genes.json \
    -s data/toxin_coding_genes_schema.json \
    -b "${SOURCE_DIR}/cdiff_blastn/cdiffamr-${genome}_blastn.tsv" \
    -o "${OUTPUT_DIR}/${genome}_toxin_coding_genes_report.json"

    #relatedness
    python3 bin/compareProfilesSingleExclude.py -f "${CGMLST_DIR}" \
    -i "${SOURCE_DIR}/cgmlst/${genome}_cgmlst.json" \
    -d 20 \
    -n 3 \
    -o "${OUTPUT_DIR}/${genome}_cgmlst_comparisons.tsv"
    
    #final report
    python3 bin/generate_individual_report.py -s "${OUTPUT_DIR}/${genome}_name.tsv" \
    -q "${OUTPUT_DIR}/${genome}_QC_summary_table.tsv" \
    -b "${OUTPUT_DIR}/${genome}_bracken_summary_table.tsv" \
    -a "${OUTPUT_DIR}/${genome}_resistance_report.json" \
    -t "${OUTPUT_DIR}/${genome}_toxin_coding_genes_report.json" \
    -r "${OUTPUT_DIR}/${genome}_cgmlst_comparisons.tsv" \
    -o "${OUTPUT_DIR}/${genome}_out_report.pdf"
done
