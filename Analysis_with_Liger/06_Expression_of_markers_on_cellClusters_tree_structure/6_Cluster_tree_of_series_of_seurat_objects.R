# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)


# Load the known cell type markers file
known_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Known_markers/CH/Hirsuta_known_cell_type_markers_latest.csv")
known_markers <- known_markers[!duplicated(known_markers$CH_ID), ]
rownames(known_markers) <- known_markers$CH_ID

# Storing directory
storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient"

clus_tree_series_of_seurat_objects(seurat_object_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_objects", 
                                   file_name_pattern = "seurat_object_of_", 
                                   marker_file = known_markers, 
                                   query_pattern = "RNA_snn_res.", 
                                   gene_ID_column = "CH_ID", 
                                   gene_name_column = "Gene_Name", 
                                   store_dir = storing_dir)

