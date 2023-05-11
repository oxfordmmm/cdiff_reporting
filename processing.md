# Summary

Uses:
* https://github.com/oxfordmmm/Bugflow_DSL2
* https://github.com/oxfordmmm/cdiff_reporting (this repo)

# Sequencing run processing

Each sequencing run needs processing woth the following steps before it can be included in the summary/ comparision reports.

## Pipleine processing

cgMLST calulcation and mapping.

```
mkdir /mnt/scratch_2/runs/test_processing/230125_VL00165_24_AACHJTHM5/cgmlst /mnt/scratch_2/runs/test_processing/230125_VL00165_24_AACHJTHM5/mapping

cd /mnt/scratch_2/runs/test_processing/230125_VL00165_24_AACHJTHM5/cgmlst
nextflow run /mnt/scratch/workflows/Bugflow_DSL2 -entry cdiff_hcgmlst_amrg_blastn_single --reads '/mnt/arun_in_bucket/UKHSA_runs/230125_VL00165_24_AACHJTHM5/*{1,2}_001.fastq.gz' --outdir /mnt/scratch_2/output/test_processing/230125_VL00165_24_AACHJTHM5/cgmlst -profile docker,oci -resume

cd /mnt/scratch_2/runs/test_processing/230125_VL00165_24_AACHJTHM5/mapping
nextflow run /mnt/scratch/workflows/Bugflow_DSL2 -entry cdiff_mapping_snpCalling_DE --reads '/mnt/arun_in_bucket/UKHSA_runs/230125_VL00165_24_AACHJTHM5/*R{1,2}_001.fastq.gz' --refFasta '/mnt/arun_in_bucket/REFSEQ_Cdiff/Cdiff_630_GCA_000009205.1.fasta' --outdir '/mnt/scratch_2/output/test_processing/230125_VL00165_24_AACHJTHM5/mapping' -profile docker,oci
```

## Pooling samples

Need to put results in general pool

* When to put in to do relatedness
* Dealing with neg samples??

```
ln -s /mnt/scratch_2/output/test_processing/230124_VL00165_23_AACHJT5M5/cgmlst/cgmlst/*.json /mnt/scratch_2/cdiff_pools_new/cgmlst_profiles/
ln -s /mnt/scratch_2/output/test_processing/230124_VL00165_23_AACHJT5M5/mapping/consensus_fa/*.fa.gz /mnt/scratch_2/cdiff_pools_new/consensus_fa/
```

## Indivdual reports

Process QC stats etc and PDF

```
cd /mnt/scratch/soft/cdiff_reporting
source venv/bin/activate
bash batch_process_cdiff.sh -s /mnt/scratch_2/output/test_processing/230125_VL00165_24_AACHJTHM5/cgmlst/ -d /mnt/scratch/soft/cdiff_reporting/data/ -c /mnt/scratch_2/cdiff_pools_new/cgmlst_profiles -o /mnt/scratch_2/output/test_processing/230125_VL00165_24_AACHJTHM5/reports
```

# Summary report batch

## Make sequence list

```
mkdir -p /mnt/scratch_2/runs/cluster_batches/230504_batch && cd /mnt/scratch_2/runs/cluster_batches/230504_batch
cp /mnt/scratch/soft/cdiff_reporting/bin/clustering_rlc/populate_seqlist.pl .
```
change to /mnt/scratch_2/runs/cluster_batches/230504_batch/consensus_fa1/. with the ls being /mnt/scratch_2/cdiff_pools_new/consensus_fa/*.gz
```
perl populate_seqlist.pl > consensus_fasta_list2.txt
```

## Make config file for runListCompare

```
cp /mnt/scratch_2/runs/cluster_batches/20230406_batch/samples.ini /mnt/scratch_2/runs/cluster_batches/230504_batch
```

change seqlist to /mnt/scratch_2/runs/cluster_batches/230504_batch/consensus_fasta_list2.txt and output_stem to /mnt/scratch_2/runs/cluster_batches/230504_batch nprocs to 10.

## Run processing

```
cd /mnt/scratch/soft/cdiff_reporting/bin/clustering_rlc
conda activate rlc_env
python runListCompare.py /mnt/scratch_2/runs/cluster_batches/230504_batch/samples.ini
```

N.B. Takes a few hours!!

## Draw cluster tress

Change dir_path of draw_all_cluster_trees.py
```
source /mnt/scratch/soft/cdiff_reporting/venv/bin/activate
cd /mnt/scratch/soft/cdiff_reporting/bin/clustering_rlc/
python3 draw_all_cluster_trees.py
```

## Make Summary PDF

```
source /mnt/scratch/soft/cdiff_reporting/venv/bin/activate
cd /mnt/scratch_2/runs/cluster_batches/230504_batch
python3 /mnt/scratch/soft/cdiff_reporting/bin/make_summary_json.py -s consensus_fasta_list2.txt -c /mnt/scratch_2/runs/cluster_batches/230504_batch -d /mnt/scratch_2/output/test_processing/230124_VL00165_23_AACHJT5M5/report /mnt/scratch_2/output/test_processing/230125_VL00165_24_AACHJTHM5/report
mkdir /mnt/scratch_2/runs/cluster_batches/230503_batch/reports
cd /mnt/scratch/soft/cdiff_reporting/bin/clustering_rlc
python3 bin/generate_summary_report.py -s /mnt/scratch_2/runs/cluster_batches/230504_batch/summary_data.json -o /mnt/scratch_2/runs/cluster_batches/230504_batch/reports
```
