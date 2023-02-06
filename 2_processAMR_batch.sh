#! /bin/bash
​
samples=("ERR3274281" "ERR3274879" "ERR3299634" "ERR3274460" "ERR3299633")
for sample in ${samples[@]}; do
    python3 2_process_AMR.py -c AMR_catalogue.json -s AMR_catalogue_schema.json -b /Users/jez/dev/cdiff_reporting/sample_data/blast/cdiffamr-${sample}.tsv -f /Users/jez/dev/cdiff_reporting/sample_data/amr_finder_plus/${sample}_forpointmuts.tsv -o /Users/jez/dev/cdiff_reporting/${sample}_resistance.json
done;