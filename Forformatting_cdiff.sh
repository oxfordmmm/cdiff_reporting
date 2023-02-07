rm QC_summary_table.tsv

#sample details eg name, mlst
perl 0_get_sample_name_MLST.pl /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/mlst/22M7591168_Leeds_S96_R_ST.tsv

#qc
perl 1_parseQCoutputs.pl /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/raw_fastqc_single/22M7591168_Leeds_S96_R_raw_reads_fastqc/22M7591168_Leeds_S96_R.txt /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/quast/22M7591168_Leeds_S96_R_Quastreport.tsv 

#amr profile
python3 AMR_process.py -c data/AMR_catalogue.json -s data/AMR_catalogue_schema.json -b /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/cdiff_blastn/cdiffamr-22M7591168_Leeds_S96_R_blastn.tsv -f /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/amr_cdiff_pointmuts/22M7591168_Leeds_S96_R_forpointmuts.tsv  -o 22M7591168_Leeds_S96_resistance_report.json

#relatedness
python3 3_compareProfilesSingleExclude.py -f /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/cgmlst/ -i /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/cgmlst/22M7591168_Leeds_S96_R_cgmlst.json -d 100 -o 22M7591168_Leeds_S96_cgmlst_comparisons.tsv

#final report
python3 AMR_report.py -s sample_name.tsv -q QC_summary_table.tsv \
-a 22M7591168_Leeds_S96_R_resistance_report.json \
-r 22M7591168_Leeds_S96_R_cgmlst_comparisons.tsv  \
-o 22M7591168_Leeds_S96_out_report.pdf

