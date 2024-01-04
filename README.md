# Clostridioides difficile Sequence Analysis Report Generator
A set of python, perl and bash scripts for generating PDF reports (of C. difficile samples) produced from the outputs of Bugflow DSL2 pipeline

## Overview
The script `batch_process_cdiff.sh` will generate all the individual sample reports. It does this by gathering the various outputs from bugflow, notably the mlst, amr/toxin profile, kraken2/bracken, cgmlst. The cgmlst is used to find related samples in the batch which is shown at the end of the report.

A python script is also run to summarise the outputs from everything into a wide summary csv.

Then to create the group report we form a list of all the consensus fastas. 
We then cluster into groups based on cgmlst distance.
For each group we use the RunListCompare process to create a phylogeny which is ten visualised with ggtree. 
Finally the group report is created giving qc pass/fail for each sample, and the cluster trees below.

## Requirements
* Python >= 3.8


## Installation
Clone the repository:
```
git clone https://github.com/oxfordmmm/cdiff_reporting
```

## Required virtual environment

First-time users:
```
cd cdiff_reporting/
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
When the virtual environmet has been populated before:
```
source  venv/bin/activate
```

You'll also need to make an environment for run list compare. The environment file is `bin/clustering_rlc/rlc_env.yaml`

Alternatively a conda environment can be created from `environment.yml`

## Usage
```
bash batch_process_cdiff.sh -s /path/to/Bugflow/output/directory/ -d ~/data \
-c /path/to/Bugflow/output/directory/cgmlst \
-o /path/to/Bugflow/output/directory/report
```

See [example_run.sh](example_run.sh) for a full example of running the BugFlow workflows and all the steps to generate all the reports.