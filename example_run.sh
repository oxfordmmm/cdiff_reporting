source /mnt/scratch/miniconda3/etc/profile.d/conda.sh

# These should match the values used for Bugflow run
batch=batch_name
run_dir=runs/$batch
output_dir=output/$batch
echo running $batch

# This should point to reporting repo
projectDir='/mnt/scratch/colpus/cdiff_reporting'

### Individual reports

# Activate reporting conda environment
conda activate $projectDir/conda_env
cd $projectDir

rm -r $output_dir/reports
bash batch_process_cdiff.sh -s $output_dir/cgmlst/ -m $output_dir/mapping -d $projectDir/data/ -c $output_dir/cgmlst/cgmlst -o $output_dir/reports

python3 $projectDir/bin/summarise_qc.py -d $output_dir/reports -o $output_dir/qc_summary.csv
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
python3 $projectDir/bin/make_cgmlst_clusters.py \
    -f $run_dir/consensus_fasta_list.tsv \
    -s $output_dir/cgmlst/cgmlst \
    -o $output_dir/cgmlst_clusters \
    -c 20

# runListCompare for each cluster (take 5 min)
cd $projectDir/bin/clustering_rlc
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
conda activate $projectDir/conda_env
Rscript $projectDir/bin/draw_trees.R $output_dir/cgmlst_clusters

# Make summary pdf
cd $run_dir
# -d can be list of directories
python3 $projectDir/bin/make_summary_json.py \
    -s consensus_fasta_list.tsv -c $output_dir/cgmlst_clusters \
    -d $output_dir/reports -o $run_dir/cluster_summary.json

cd $projectDir
python3 bin/generate_summary_report.py -s $run_dir/cluster_summary.json -o $output_dir/cluster_report.pdf

