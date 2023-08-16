# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Convert the liger object to a seurat object to perform downstream analysis
###

storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient"

liger_to_seurat_conversion(liger_object_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/Liger_objects", 
                           file_name_pattern = "Liger_object_K_", 
                           scale_coefficient = FALSE, 
                           store_dir = storing_dir,
                           store_folder = "Seurat_objects")
