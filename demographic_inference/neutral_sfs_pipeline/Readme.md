# Pipeline to build the neutral SFS

## 1. Neutral_PreviewProjection_EasySFS_1D.sh
Identifies variants and invariants sites for each chromosome/scaffold vcf file. Additionally, it generates the projection of the folded SFS using the easySFS program and plots the projection for each chromosome/scaffold. Calls and uses other scripts and files: easySFS_a.py, easySFS_function_determineOptimalProjection.R and pop_map.txt
  
## 2. SFS_preview_v1.R
Puts together the preview of SFS projections for each chromosomes/scaffold

## 3. SFS_projection_chr.sh
Runs the SFS projection per chromosomes. Calls and uses other scripts and files: easySFS_a.py, pop_map.txt
  
## 4. Projection_file_list.sh
This script concentrate files names of SFS projection per chromosomes

## 5. SFS_projection_visualization.R
Join SFS projections for all the chromosomes and plot them

## 6. SFS_count_monomorphic_sites.sh (calls other scripts and files)
  getMonomorphicProjectionCounts.1D.2DSFS.py
  
  pop_map.txt
  
## 7. SFS_projection_visualization_mono.R


