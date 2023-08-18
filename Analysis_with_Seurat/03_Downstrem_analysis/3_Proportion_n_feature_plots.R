# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Function to generate results for a seurat object with clustering results for multiple resolution parameter values
resolution_explorer <- function(seuratObject, 
                                store_dir = NULL, 
                                store_folder = "Your_choice_of_name"
){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # let's go over the resolution values
  temp = str_sort(names(seuratObject@meta.data)[str_detect(names(seuratObject@meta.data), pattern = "res")], numeric = TRUE)
  
  for (i in c(1:length(temp))){
    
    temp_str = substr(gsub("[^[:alnum:]]", "", temp[i]), start = nchar(gsub("[^[:alnum:]]", "", temp[i])) - 1, stop = nchar(gsub("[^[:alnum:]]", "", temp[i])))
    temp_str = str_c("RES_", gsub(pattern = "[A-Za-z]", replacement = "", temp_str))
    
    if (!dir.exists(str_c(temp_dir, "/", temp_str))){
      dir.create(str_c(temp_dir, "/", temp_str), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_res_dir = str_c(temp_dir, "/", temp_str, "/")
    
    # Set the cluster identity to the seurat object
    Idents(seuratObject) <- temp[i]
    
    Idents(seuratObject) <- factor(Idents(seuratObject), levels = seq(0, length(levels(Idents(seuratObject))) - 1))
    
    replicate_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotype")
    
    species_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotype", rep_prop = FALSE)
    
    species_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotype", rep_prop = TRUE)
    
    feature_count_n_proportion(seuratObject = seuratObject, store_dir = temp_res_dir, gene_ID = "Chir02Ox-b03100.2", gene_name = "STM")
    
    feature_count_n_proportion(seuratObject = seuratObject, store_dir = temp_res_dir, gene_ID = "Chir06Ox-b35150.2", gene_name = "RCO")
    
    RDimension_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Reduced_representation",dimension_reduction_name = "umap")
    
    # RDimension_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Reduced_representation",dimension_reduction_name = "tsne")
    
    genes_feature_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Genes_feature_plot", genes_ID = c("Chir02Ox-b03100.2", "Chir06Ox-b35150.2"), genes_name = c("STM", "RCO"), reduction_name = "umap")
    
    genes_feature_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Genes_feature_plot", genes_ID = c("Chir02Ox-b03100.2", "Chir06Ox-b35150.2"), genes_name = c("STM", "RCO"), reduction_name = "tsne")
    
    markers_dot_plot(seuratObject = seuratObject, store_dir = temp_res_dir, marker_file = cell_type_markers, genes_ID_column = 4, genes_name_column = 1, group_clusters = TRUE)
    
    # side_by_side_proportion_comparison(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable_name = "Species")
    
    RDimension_split_group_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Species")
    
    # cells_per_cluster_split_viz(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Species")
    
    feature_expression_split_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Species", gene_ID = "Chir06Ox-b35150.2", gene_name = "RCO")
    
    feature_expression_split_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Species", gene_ID = "Chir02Ox-b03100.2", gene_name = "STM")
    
    # cluster_conserved_marker_finder(seuratObject = seuratObject, store_dir = temp_res_dir, grouping_variable = "Species")
    #
    # between_groups_differentially_expressed_genes(seuratObject = seuratObject, store_dir = temp_res_dir, between_group_variable = "Species")
    # 
    # group_n_cluster_DEG_finder(seuratObject = seuratObject, store_dir = temp_res_dir)
    # 
    # cluster_marker_finder(seuratObject = seuratObject, store_dir = temp_res_dir)
    
  }
  
}

# Load
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_rco_seurat.RData")

# Storing directory
storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1"

cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Curated_markers_hirsuta.csv")

resolution_explorer(seuratObject = integrated.data, store_dir = storing_dir, store_folder = "Results_for_a_range_of_resolution")
