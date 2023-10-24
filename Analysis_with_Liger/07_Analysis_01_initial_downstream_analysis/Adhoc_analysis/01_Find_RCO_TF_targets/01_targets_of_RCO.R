project_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_object/seurat_object_of_K_50.RData")

# Set the identity of the cells to a particular resolution value
Idents(integrated.data) <- integrated.data$RNA_snn_res.0.3

# In the metadata container add the cell ID labels
integrated.data$cell_ID = rownames(integrated.data@meta.data)

# Check the RCO expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", ] != 0])

# Get the details for the RCO-expressing cells
rco_cells_details = integrated.data@meta.data[Cells_with_expression_detection, ]

# Subset the seurat object to keep only those cells with RCO's expression detected 
rco_cells_subset = subset(integrated.data, subset = cell_ID %in% Cells_with_expression_detection)

# Combine all the genes that was identified as specific markers or GEP genes for the RCO cell type in the analysis of WT C.hirsuta

# Characteristic GEP, GEP 32, of the RCO cell type
GEP_genes = read.delim(paste0(project_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_output/Factor_Ks_GEPs/Factor_K_44_GEPs/GEP_32.txt"), header = FALSE)[, 1]

# Specific markers of the RCO cell type in the Liger analysis
rco_markers_liger = read.csv(paste0(project_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Cluster_specific_markers/Cluster_5_specific_markers.csv"))
rco_markers_liger = rco_markers_liger$gene_ID

# Specific markers of the RCO cell type in the Liger analysis
rco_markers_seurat = read.csv("seurat_rco_cluster_7.csv")
rco_markers_seurat = rco_markers_seurat$gene_ID

# Combine the RCO cell type markers
combined_rco_markers = unique(c(GEP_genes, rco_markers_liger, rco_markers_seurat)) # 124 genes

writeLines(combined_rco_markers, "all_rco_markers.txt")
