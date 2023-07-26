source /mnt/scratch/miniconda3/etc/profile.d/conda.sh

batch="test_batch"

input='/mnt/arun_in_bucket/UKHSA_runs/230125_VL00165_24_AACHJTHM5/*{1,2}_001.fastq.gz'
ref='/mnt/arun_in_bucket/REFSEQ_Cdiff/Cdiff_630_GCA_000009205.1.fasta'
run_dir=/mnt/scratch_2/runs/colpus_test/$batch
output_dir=/mnt/scratch_2/output/colpus_test/$batch


mkdir -p $run_dir/cgmlst $run_dir/mapping

cd $run_dir/cgmlst
nextflow run /mnt/scratch/colpus/Bugflow_DSL2 -entry cdiff_hcgmlst_amrg_blastn_single --reads $input --outdir $output_dir/cgmlst -profile docker,oci -resume

cd $run_dir/mapping
nextflow run /mnt/scratch/colpus/Bugflow_DSL2 -entry cdiff_mapping_snpCalling_DE --reads $input --refFasta $ref --outdir $output_dir/mapping -profile docker,oci

### Individual reports
conda activate /mnt/scratch/colpus/cdiff_reporting/conda_env
cd /mnt/scratch/colpus/cdiff_reporting

rm -r $output_dir/reports
bash batch_process_cdiff.sh -s $output_dir/cgmlst/ -m $output_dir/mapping -d /mnt/scratch/colpus/cdiff_reporting/data/ -c $output_dir/cgmlst/cgmlst -o $output_dir/reports

python3 /mnt/scratch/colpus/cdiff_reporting/bin/summarise_qc.py -d $output_dir/reports -o $output_dir/qc_summary.csv
### For group reports

# Make consensus fasta list
> $run_dir/consensus_fasta_list.tsv
for fasta in $output_dir/mapping/consensus_fa/*
do
    basename=$(basename $fasta)
    simplename=${basename/.fa.gz/}

    if [[ $simplename != neg* ]]; then
        echo -e "$simplename\t$fasta" >> $run_dir/consensus_fasta_list.tsv
    fi
done


# Make cgmlst clusters
mkdir $output_dir/cgmlst_clusters 
python3 /mnt/scratch/colpus/cdiff_reporting/bin/make_cgmlst_clusters.py \
    -f $run_dir/consensus_fasta_list.tsv \
    -s $output_dir/cgmlst/cgmlst \
    -o $output_dir/cgmlst_clusters \
    -c 20

# runListCompare for each cluster (take 5 min)
cd /mnt/scratch/colpus/cdiff_reporting/bin/clustering_rlc
conda activate rlc_env

for cluster in $output_dir/cgmlst_clusters/*consensus.tsv
do
    base=$(basename $cluster)
    i=$(echo $base | cut -d'_' -f2)
    echo cluster $i
    echo $cluster
    python runListCompare.py -i runListCompare_config_base.ini -s $cluster -o $output_dir/cgmlst_clusters/$i 

    for tree in $output_dir/cgmlst_clusters/$i/cluster_ml/*scaled*.tree
    do
        echo $tree
        type=$(basename $tree | cut -d'_' -f3)
        cp $tree $output_dir/cgmlst_clusters/cluster_${i}_${type}_scaled.newick
    done
done

# Draw all trees
conda activate /mnt/scratch/colpus/cdiff_reporting/conda_env
Rscript /mnt/scratch/colpus/cdiff_reporting/bin/draw_trees.R $output_dir/cgmlst_clusters

# Make summary pdf
cd $run_dir
# -d can be list of directories
python3 /mnt/scratch/colpus/cdiff_reporting/bin/make_summary_json.py \
    -s consensus_fasta_list.tsv -c $output_dir/cgmlst_clusters \
    -d $output_dir/reports -o $run_dir/cluster_summary.json

cd /mnt/scratch/colpus/cdiff_reporting
python3 bin/generate_summary_report.py -s $run_dir/cluster_summary.json -o $output_dir/cluster_report.pdf

