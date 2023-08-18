# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)


# Load the known cell type markers file
known_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Hirsuta_known_cell_type_markers_latest.csv")

# Load
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_rco_seurat.RData")

DefaultAssay(integrated.data) <- "RNA"

# Storing directory
storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1"

clus_tree_of_cell_clusters(seuratObject = integrated.data, 
                           marker_file = known_markers, 
                           query_pattern = "integrated_snn_res.", 
                           gene_ID_column = "CH_ID", 
                           gene_name_column = "Gene_Name", 
                           store_dir = storing_dir)
