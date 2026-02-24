#!/bin/bash
# Tree input generation for scRNA-seq data

set -e

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
echo " ----- Start at : $startTime -----"

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Load configuration
CONFIG_FILE="${PROJECT_ROOT}/config/paths.yaml"
if [ -f "$CONFIG_FILE" ]; then
    RNA_EDITING_FILE=$(grep 'rna_editing_file:' "$CONFIG_FILE" | grep -v '#' | head -1 | awk '{print $2}' | tr -d '"')
else
    RNA_EDITING_FILE=${RNA_EDITING_FILE:-""}
fi

# Add script directory to Python path
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH}"

# Parse command line arguments
export workpath=$1
cd ${workpath}
export sample1=$2
export barcode1=$3
export bam_file1=$4
export mutation_list=$5
export cellnum=$6
export thread=$7

echo "===== Tree input generation parameters ====="
echo "workpath: ${workpath}"
echo "sample1: ${sample1}"
echo "barcode1: ${barcode1}"
echo "bam_file1: ${bam_file1}"
echo "mutation_list: ${mutation_list}"
echo "cellnum: ${cellnum}"
echo "thread: ${thread}"
echo "============================================="

echo "===== Start running: ====="

################################################################
# Generate tree input raw files
echo "### Step1: Generate tree input raw files ###"
mkdir -p ${workpath}/treeinput/

python -m others.prepare_phylo_input_from_scRNA \
    --barcode_files ${barcode1} \
    --bams ${bam_file1} \
    --mutlist ${mutation_list} \
    --outprefix ${workpath}/treeinput/treeinput \
    --samples ${sample1}

################################################################
# Data preprocess before phylogeny
echo "### Step2: Data preprocess before phylogeny ###"
Rscript ${SCRIPT_DIR}/1.rawPreprocess_spatial.extract_identifier_sites.R \
    --inputfile ${workpath}/treeinput/treeinput_spot_c_${cellnum}.csv \
    --cellnum ${cellnum} \
    --outputpath ${workpath}/data \
    --scid_file ${workpath}/treeinput/treeinput_scid_barcode.txt \
    --is_remove_cells no \
    --threshold 0.5 \
    --indid ${sample1}

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
echo " ----- End at : $endTime_s -----"
sumTime=$[ $endTime_s - $startTime_s ]
echo "$startTime ---> $endTime" "Total:$sumTime seconds"