ls /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/assemblies/*_R_contigs.fa > processed_genomes.tsv
perl -p -i -e 's/\_contigs\.fa//g' processed_genomes.tsv
perl -p -i -e 's/\/mnt\/scratch\_2\/output\/out1\_230125\_VL00165\_24\_AACHJTHM5\/assemblies\///g' processed_genomes.tsv

for genome in $(cat processed_genomes.tsv)
do
rm QC_summary_table.tsv
#sample details eg name, mlst
perl 0_get_sample_name_MLST.pl /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/mlst/$genome\_ST.tsv
#qc
perl 1_parseQCoutputs.pl /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/raw_fastqc_single/$genome\_raw\_reads\_fastqc/$genome.txt /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/quast/$genome\_Quastreport.tsv 
#amr profile
python3 2_process_AMR_profiles.py -c data/AMR_catalogue.json \
-s data/AMR_catalogue_schema.json \
-b /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/cdiff_blastn/cdiffamr-$genome\_blastn.tsv \
-f /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/amr_cdiff_pointmuts/$genome\_forpointmuts.tsv  \
-o $genome\_resistance_report.json
#relatedness
python3 3_compareProfilesSingleExclude.py -f /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/cgmlst/ \
-i /mnt/scratch_2/output/out1_230125_VL00165_24_AACHJTHM5/cgmlst/$genome\_cgmlst.json \
-d 50 \
-o $genome\_cgmlst_comparisons.tsv
#final report
#python3 4_create_Cdiff_report.py -s sample_name.tsv \
#-q QC_summary_table.tsv -a $genome\_resistance_report.json \
#-r $genome\_cgmlst_comparisons.tsv \
#-o $genome\_out_report.pdf
done