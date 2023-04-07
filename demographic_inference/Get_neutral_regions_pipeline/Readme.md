# Pipeline to get neutral regions to build neutral SFS

## 1. GetHiQualCoords_20200602.sh

Obtains the coordinates of high quality sites from the final VCF file out of the last step of the variant calling pipeline (WGSproc8). This bash script runs the script obtain_high_qual_coordinates.py that only uses sites that passed all filters of our variant calling pipeline (i.e. have PASS in the FILTER column), but also eliminates sites that do not have information in the INFO column (i.e. it shows only a dot "."). This happens because GATK variant filtration will annotate a site as PASS when there is no information in the INFO column. It uses the script "obtain_high_qual_coordinates.py"


## 2. exon_distance_calc.sh

It uses the file "all_HQCoords_sorted_merged.bed" resulting from the previous script to calculate different distances from exons and coding regions, the results from this two features might be different or the same depending on the quality of your annotation. The results of this script will help you to decide what distance to use to define the neutral regions.


## 3. Extract_noCpG_noRepetitive.sh

Extracts only the high quality genomic regions that are at least at certained predefined distance from exons or coding regions (in our case was 20Kb to define our neutral regions) and not in CpG islands or repetitive regions. To do this it uses the files that were generated with the concatenation procedure of repepetitive and cpg enriched regions obtained after running the "concat_bed_files_steps.txt" file and the "all_HQCoords_min20kb_DistFromCDR.0based.bed" file generated after running the "exon_distance_calc.sh" script. This filter files are applied on the final VCF files from last step of the variant calling pipeline (WGSproc8) "JointCalls_all50_08_B_VariantFiltration_${IDX}.vcf.gz"

## 4. get_Coord_file.sh 

Gets the neutral coordinates without cpg or repetitive regions from the vcf files generated in the previous step ("nocpg_repeats_SFS_${IDX}.vcf.gz") to bed files. It uses the script "obtain_noCpG_noRepetitive_coordinates.py"

## 5. identify_Conserved_Regions.sh

Identifies regions that match to conserved regions (aligns to zebra fish genome), processing through the befiles generated in the previous step. The conserved regions are removed from further analysis in the next step.

## 6. neutral_Sites2vcf.sh

Extracts the neutral regions that are not conserved (i.e. do not map to zebra fish genome), to new vcf files per chromosome or scaffold. These extracted neutral regions will be used to buil the netural SFS projection. 


