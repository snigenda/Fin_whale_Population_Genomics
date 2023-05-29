# Pipeline to build the neutral SFS

## 1. Neutral_PreviewProjection_EasySFS_1D.sh (calls other scripts and files)
This scripts identifies variants and invariants sites for specified vcf files. Additionally, generates and plots the folded SFS preview using the sites identified. This scripts calls and uses other scripts and files: easySFS_a.py, easySFS_function_determineOptimalProjection.R and pop_map.txt
  
## 2. SFS_preview_v1.R
## 3. SFS_projection_chr.sh (calls other scripts and files)
  easySFS_a.py
  
  pop_map.txt
  
## 4. Projection_file_list.sh
## 5. SFS_projection_visualization.R
## 6. SFS_count_monomorphic_sites.sh (calls other scripts and files)
  getMonomorphicProjectionCounts.1D.2DSFS.py
  
  pop_map.txt
  
## 7. SFS_projection_visualization_mono.R


