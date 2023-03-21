# Clostridioides difficile Sequence Analysis Report Generator
A set of python, perl and bash scripts for generating PDF reports (of C. difficile samples) produced from the outputs of Bugflow DSL2 pipeline

## Requirements
* Python >= 3.8


## Installation
Clone the repository:
```
git clone https://github.com/oxfordmmm/cdiff_reporting
```

Required virtual environment:
```
cd cdiff_reporting/
source  venv/bin/activate
```

## Usage
```
bash batch_process_cdiff.sh -s /path/to/Bugflow/output/directory/ -d ~/data \
-c /path/to/Bugflow/output/directory/cgmlst \
-o /path/to/Bugflow/output/directory/report
```