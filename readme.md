# Whole Genome Sequence and Assemly

## Overview
## Workflow
### 1. Raw Reads

 We will use grabseqs for raw reads downloading.
#### 1.1. Install grabseqs

```bash
# install grabseqs
conda create -n grabseqs -y
conda activate grabseqs
conda install python=3.9 -y
pip install grabseqs
# dependencies
conda install conda-forge::pigz -y
conda install bioconda::sra-tools -y
```

#### 1.2. Download raw reads

```bash
# to download a sequence illumina run
grabseqs sra -t 4 -m metadata.csv SRR8893090

# first run for nanopor reads
grabseqs sra -t 4 -m metadata.csv SRR8893087
# second run for nanopore reads
grabseqs sra -t 4 -m metadata.csv SRR8893086
# pacbio reads
grabseqs sra -t 4 -m metadata.csv SRR8893091
```
---

