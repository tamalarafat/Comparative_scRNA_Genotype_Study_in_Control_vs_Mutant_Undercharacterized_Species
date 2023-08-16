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
    
    replicate_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotypes")
    
    group_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotypes", rep_prop = FALSE)
    
    group_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotypes", rep_prop = TRUE)
    
    feature_count_n_proportion(seuratObject = seuratObject, store_dir = temp_res_dir, gene_ID = "Chir02Ox-b03100.2", gene_name = "STM")
    
    feature_count_n_proportion(seuratObject = seuratObject, store_dir = temp_res_dir, gene_ID = "Chir06Ox-b35150.2", gene_name = "RCO")
    
    RDimension_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Reduced_representation",dimension_reduction_name = "umap")
    
    genes_feature_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Genes_feature_plot", gene_IDs = c("Chir02Ox-b03100.2", "Chir06Ox-b35150.2"), genes_name = c("STM", "RCO"), reduction_name = "umap")
    
    genes_feature_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Genes_feature_plot", gene_IDs = c("Chir02Ox-b03100.2", "Chir06Ox-b35150.2"), genes_name = c("STM", "RCO"), reduction_name = "tsne")
    
    markers_dot_plot(seuratObject = seuratObject, store_dir = temp_res_dir, marker_file = cell_type_markers, genes_ID_column = 4, genes_name_column = 2, group_clusters = TRUE)
    
    RDimension_split_group_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes")
    
    cells_per_cluster_split_viz(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes")
    
    feature_expression_split_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes", gene_ID = "Chir06Ox-b35150.2", gene_name = "RCO")
    
    feature_expression_split_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes", gene_ID = "Chir02Ox-b03100.2", gene_name = "STM")
    
  }
  
}

# For a series of seurat objects
resolution_explorer_series_of_seurat_objects <- function(seurat_object_dir, 
                                                         file_name_pattern, 
                                                         store_dir = NULL,
                                                         store_folder = "Cells_genes_clusters",
                                                         ident_levels = NULL,
                                                         add_replicate_labels = NULL,
                                                         add_species_labels = NULL,
                                                         add_tissue_labels = NULL
){
  
  # Lets get the saved file names
  seuratObjects = str_sort(list.files(path = seurat_object_dir, pattern = file_name_pattern), numeric = TRUE)
  
  # Let's create a directory to store the GO annotation results for different factorization
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Let's create a directory to store the GO annotation results
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_main_dir = str_c(store_dir, "/", store_folder, "/")
  
  for (i in c(1:length(seuratObjects))){
    # Creating the file name with path location
    file_name = str_c(seurat_object_dir, "/", seuratObjects[i])
    
    # Loading the liger object and storing it to a temporary variable
    integrated.data = loadRData(file_name) # Load the RData file and assign it to a variable using the function loadRData
    
    # Dataset information
    integrated.data$Replicates <- integrated.data$orig.ident
    integrated.data$Replicates <- factor(integrated.data$Replicates, levels = c("M3", "M6", "M7", "O1", "O2", "O3", "O7"), labels = c("rco-3", "rco-6", "rco-7", "OX-1", "OX-2", "OX-3", "OX-7"))
    
    # Genotypes information
    integrated.data$Genotypes <- integrated.data$orig.ident
    integrated.data$Genotypes <- substr(integrated.data$Genotypes, 1, nchar(as.character(integrated.data$Genotypes)) - 1)
    integrated.data$Genotypes <- factor(integrated.data$Genotypes, levels = c("M", "O"), labels = c("mutant", "WT"))
    
    # Tissue information
    integrated.data$Tissue <- integrated.data$orig.ident
    integrated.data$Tissue <- substr(integrated.data$Tissue, 1, nchar(as.character(integrated.data$Tissue)) - 1)
    integrated.data$Tissue <- factor(integrated.data$Tissue, levels = c("M", "O"), labels = c("Leaf", "Leaf"))
    
    # Save the object with the added metadata information
    save(integrated.data, file = str_c(seurat_object_dir, "/", "seurat_object_of_K_", dim(integrated.data@reductions[["inmf"]])[2], ".RData"))
    
    resolution_explorer(seuratObject = integrated.data, store_dir = temp_main_dir, store_folder = str_c("Factor_K_", ncol(integrated.data@reductions$inmf)))
    
  }
}

# Known cell type markers
cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Known_markers/CH/Curated_markers_hirsuta.csv")

# Directory containing the seurat objects
seurat_files_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_objects"

# Storing directory
storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient"

resolution_explorer_series_of_seurat_objects(seurat_object_dir = seurat_files_dir, 
                                             file_name_pattern = "seurat_object_of_", 
                                             store_dir = storing_dir, 
                                             store_folder = "Explore_resolutions")



