project_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = paste0(project_dir, "comparative_study_of_wt_and_mutant_genotype_Cardamine/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs/Differentially_expressed_genes/DEGs_between_genotypes_of_a_cluster/DEG_between_Genotypes_DEtest_wilcox")

# Storing directory
storing_dir = paste0(project_dir, "comparative_study_of_wt_and_mutant_genotype_Cardamine/Analysis_with_Liger/07_Analysis_01_initial_downstream_analysis/Analysis_outputs")

specific_markers_list = specific_marker_finder(DEG_file_dir = DEG_dir, pct_detection = 0.1, pct_diff = 0.3, store_dir = storing_dir)