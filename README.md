# PhyloSOLID

Tree building from single-cell sequencing data (scRNA-seq and scDNA-seq)

## Overview

PhyloSOLID is a comprehensive pipeline for building phylogenetic trees from single-cell sequencing data. It supports both scRNA-seq and scDNA-seq modes, with features for mutation filtering, tree construction, and visualization.

## Installation

### Prerequisites

- Python 3.7 or higher
- R 4.0 or higher with required packages
- samtools, bedtools
- ANNOVAR (for variant annotation)

### Step 1: Clone the repository

```bash
git clone https://github.com/yourusername/phylosolid.git
cd phylosolid
```

### Step 2: Install Python dependencies
```bash
pip install -r requirements.txt
```

### Step 3: Install in development mode
```bash
pip install -e .
```

### Step 4: Install ANNOVAR

##### PhyloSOLID uses ANNOVAR for variant annotation. You need to install it separately:

1. Register and download ANNOVAR:
	Visit ANNOVAR Download Page
	Fill in the registration form (free for academic use)
	Download the latest version (annovar.latest.tar.gz)

2. Install ANNOVAR:
```
tar -xzvf annovar.latest.tar.gz -C /path/to/software/
cd /path/to/software/annovar
```

3. Download required databases (hg38 build):
```
# Create humandb directory
mkdir -p /path/to/software/annovar/humandb

# Download refGene database (required)
./annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/

# Download additional databases (recommended)
./annotate_variation.pl -downdb -buildver hg38 -webfrom annovar cytoBand humandb/
./annotate_variation.pl -downdb -buildver hg38 -webfrom annovar snp138 humandb/
```

### Step 5: Configure paths
##### Copy the template configuration file and update it with your paths:
```
cp config/paths.yaml.template config/paths.yaml
```
##### Edit config/paths.yaml to set the correct paths for:
	ANNOVAR installation directory
	Reference files (genome FASTA, GFF3, etc.)
	Other external tools

### Step 6: Verify installation
```
# Check if ANNOVAR is properly configured
phylosolid check-annovar --config config/paths.yaml

# Run the test pipeline
cd demo
./run_demo.sh
```

## Usage
### scRNA-seq mode
##### Basic usage:
```
phylosolid scrna \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt
```
##### With all options:
```
phylosolid scrna \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt \
    --metadata metadata.txt \
    --read-len 91 \
    --cellnum 155 \
    --threads 8 \
    --running-type production \
    --workdir ./results \
    --config config/paths.yaml
```

### Running specific steps
```
# Run only feature extraction
phylosolid scrna --sample SAMPLE_ID ... --steps feature_extraction

# Run only tree input generation
phylosolid scrna --sample SAMPLE_ID ... --steps tree_input

# Run only tree building
phylosolid scrna --sample SAMPLE_ID ... --steps tree_building
```

### Parallel execution
##### Run feature extraction and tree input in parallel:
```
phylosolid scrna --sample SAMPLE_ID ... --parallel
```

### scDNA-seq mode
```
phylosolid scdna \
    --sample SAMPLE_ID \
    --mutation-list mutations.txt \
    --bam sample.bam \
    --barcode barcodes.txt
```

## Configuration
### Paths configuration (config/paths.yaml)
```
# ANNOVAR configuration
annovar:
  script_dir: "/path/to/annovar"      # Directory containing annotate_variation.pl
  humandb: "/path/to/annovar/humandb" # Directory with downloaded databases
  build: "hg38"                        # Genome build (hg38/hg19)

# Reference files
reference:
  genome_fasta: "/path/to/genome.fa"
  gff3_file: "/path/to/wgEncodeGencodeExonSupportV44.sort.bed"
  mappability_file: "/path/to/k24.umap.bedgraph"
  gnomad_file: "/path/to/hg38_gnomad312_genome_only_af_all.txt"
  rna_editing_file: "/path/to/COMBINED_RADAR_REDIprotal_DARNED_hg38_all_sites.bed"
  run_get_prior: "/path/to/run_get_prior.py"
```

## Input File Formats
### Mutation list (mutations.txt)
```
chr1_1000_A_G_gene1
chr1_2000_C_T_gene2
chr2_3000_G_A_gene1
chr3_4000_T_C_gene3
```
Format: chromosome_position_reference_alt_gene

### Barcode file (barcodes.txt)
```
AAACCTGAGAAACCAT-1
AAACCTGAGAAACCGG-1
AAACCTGAGAAACCTA-1
```

### Metadata file (metadata.txt)
```
cell_barcode    cell_type
AAACCTGAGAAACCAT-1  CD8+T
AAACCTGAGAAACCGG-1  CD4+T
AAACCTGAGAAACCTA-1  Monocyte
```

## Output Structure
After running, results are organized as:
```
workdir/
└── SAMPLE_ID/
    ├── 01_features/              # Feature extraction results
    │   ├── depth_in_spots/
    │   ├── SAMPLE_ID.benchmark_patched.feature.txt
    │   └── ...
    ├── 02_treeinput/              # Tree input files
    │   ├── treeinput/
    │   │   ├── treeinput_spot_c_155.csv
    │   │   ├── treeinput_scid_barcode.txt
    │   │   └── features_file.txt
    │   └── data/                   # Preprocessed data
    ├── 03_tree_building/           # Tree building results
    │   └── results/
    │       └── PhyloSOLID/
    │           ├── cell_by_mut.CFMatrix
    │           └── celltree.newick
    └── pipeline_summary.yaml       # Pipeline execution summary
```


## Demo
A complete demo with test data is available in the demo/ directory:
```
cd demo
./run_demo.sh
```

This will:

1. Create test input files
2. Run the full pipeline on test data
3. Validate the output

### Troubleshooting

### ANNOVAR not found
If you see "ANNOVAR not found" error:

1. Check that ANNOVAR is installed
2. Verify the paths in config/paths.yaml
3. Run phylosolid check-annovar to diagnose

### Read length mismatch
If you see warnings about read length:

1. Check your BAM file's actual read length:
```
samtools view your.bam | head -n1 | awk '{print length($10)}'
```
2. Use the --read-len parameter to specify the correct value

### Permission denied errors
Ensure scripts are executable:
```
chmod +x scripts/scrna/**/*.sh
chmod +x scripts/scrna/**/*.py
chmod +x scripts/scrna/**/*.R
```


## Citation

If you use PhyloSOLID in your research, please cite:

    1.Yang, Q. et al. PhyloSOLID: robust phylogeny reconstruction from single-cell data despite pervasive errors and extreme sparsity.

## License

MIT License

Copyright (c) 2026 douymLab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


