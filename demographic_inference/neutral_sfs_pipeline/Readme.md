# Pipeline to build the neutral SFS

## 1. Neutral_PreviewProjection_EasySFS_1D.sh (uses easySFS_a.py and pop_map.txt)
## 2. SFS_preview_v1.R
## 3. SFS_projection_chr.sh (uses easySFS_a.py and pop_map_filtered.txt)
## 4. Projection_file_list.sh
## 5. SFS_projection_visualization.R
## 6. SFS_count_monomorphic_sites.sh (calls getMonomorphicProjectionCounts.1D.2DSFS.py and pop_map_filtered.txt)
## 7. SFS_projection_visualization_mono.R

easySFS_1_ProjectionPreview.sh                 pop_map.txt
easySFS_ab.py                                  Projection_file_list.sh
easySFS_a.py                                   Projection_file_list.sh.save
easySFS_function_determineOptimalProjection.R  SFS_preview_v1.R
Fin_whale_pipeline                             SFS_projection_chr.sh
Neutral_PreviewProjection_EasySFS_1D.sh        SFS_projection_visualization.R
pop_map_filtered.txt

