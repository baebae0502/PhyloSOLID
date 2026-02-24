#!/bin/bash

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
echo " ----- Start at : $startTime -----"

##### load modules
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate /home/douyanmeiLab/yangqing/Conda/envs/pmg
module load samtools
module load bedtools
module load bcftools

##### Setting and loading
export workpath=$1
cd ${workpath}
export sample1=$2
export cellnum=$3
export celltype_file=$4


# workpath="/storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/samples/phylo_dna_1465/"
# cd ${workpath}
# sample1="UMB1465"


echo "===== Check the parameters ====="
echo "workpath: "${workpath}
echo "sample1: "${sample1}
echo "cellnum: "${cellnum}
echo "celltype_file: "${celltype_file}

echo "===== Start running: ====="


################################################################
# Data preprocess before phylogeny
echo "##### Data preprocess before phylogeny"
echo "${workpath}/treeinput is ready"
echo "${workpath}/data will be generated"

Rscript /storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/1.rawPreprocess_yunchao.R \
  --inputfile ${workpath}/treeinput/treeinput_file.txt \
  --cellnum ${cellnum} \
  --outputpath ${workpath}/data \
  --scid_file ${workpath}/treeinput/treeinput_scid_barcode.txt \
  --is_remove_cells no \
  --threshold 0.9 \
  --indid ${sample1}


################################################################
# 1. PhyloSOLID
echo "##### Run PhyloSOLID for tree"
echo "${workpath}/results will be generated"

cd /storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src
python -m phylosolid.run_phylosilid_fullTree_scDNA \
--sampleid ${sample1} \
--inputpath ${workpath}/data \
--outputpath ${workpath}/results \
--celltype_file ${celltype_file} \
--features_file ${workpath}/treeinput/features_file.txt \


cd ${workpath}/results
mkdir -p ${workpath}/results/PhyloSOLID
cd ${workpath}/results/PhyloSOLID
# cp ${workpath}/results/mutation_integrator/phylo/final_cleaned_M_full_basedPivots.filtered_sites_inferred.CFMatrix ${workpath}/results/PhyloSOLID/cell_by_mut.CFMatrix
cp ${workpath}/results/scaffold_builder/phylo_scaffold_tree/final_cleaned_M_scaffold_basedPivots.filtered_sites_inferred.CFMatrix ${workpath}/results/PhyloSOLID/cell_by_mut.CFMatrix

Rscript /storage/douyanmeiLab/yangqing/tools/PhyloMosaicGenie/pmg/src/phylosolid/all_methods/convert_PhyloSOLID_tree.R ${workpath}/results/PhyloSOLID/cell_by_mut.CFMatrix ${workpath}/results/PhyloSOLID/celltree.newick






endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`
echo " ----- End at : $startTime ----> $endTime -----"
sumTime=$[ $endTime_s - $startTime_s ]
echo "Total: $sumTime seconds"