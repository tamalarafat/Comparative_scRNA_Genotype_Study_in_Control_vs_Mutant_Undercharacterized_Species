# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_marker_identification", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_matrix_manipulation", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_object/seurat_object_of_K_50.RData")

# Set the identity of the cells to a particular resolution value
Idents(integrated.data) <- integrated.data$RNA_snn_res.0.3

# In the metadata container add the cell ID labels
integrated.data$cell_ID = rownames(integrated.data@meta.data)

# Check the RCO expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", ] != 0])

# Subset the seurat object
cells_subset = subset(integrated.data, subset = cell_ID %in% Cells_with_expression_detection)


# Subset the gene set
# Genes to keep - From the previous analysis of WT Cardamine, get the characteristic GEP of the RCO cell type (GEP32), and RCO cell type markers from the Liger (cluster 5) and Seurat (cluster 7) analysis

GEP_genes = read.delim("RCO_cluster_markers_from_wt_cardamine_analysis/GEP_genes/GEP_32.txt", header = FALSE)[, 1]

markers_seurat = read.csv("RCO_cluster_markers_from_wt_cardamine_analysis/Markers_seurat_analysis/Cluster_7_specific_markers.csv", row.names = 1)[, "gene_ID"]

markers_liger = read.csv("RCO_cluster_markers_from_wt_cardamine_analysis/Markers_liger_analysis/Cluster_5_specific_markers.csv")[, "gene_ID"]


# How many common genes are present in these markers lists
paste0("Among the ", length(GEP_genes), " GEP genes, ", length(intersect(GEP_genes, markers_liger)), " and ", length(intersect(GEP_genes, markers_seurat)), " genes are present in the specific markers sets identfied for the RCO cell type in the Liger and Seurat analysis, respectively.")
paste0(length(intersect(markers_seurat, markers_liger)), " common markers between the specific markers sets obtained in Liger (", length(markers_liger), ") and Seurat (", length(markers_seurat), ") analysis.")

# Combine the markers sets
combined_markers = unique(c(GEP_genes, markers_liger, markers_seurat)) # 124 genes

# Save the combined list
writeLines(combined_markers, "RCO_cluster_all_analysis_markers.txt")

# Perform differential expression analysis between the two groups or conditions (WT and rco-mutant)

# Set the idents to the genotype information so we know the two groups of cells
Idents(cells_subset) <- "Genotypes"

# Perform DE analysis
degs_of_subsets = FindMarkers(cells_subset, ident.1 = "WT", features = combined_markers, test.use = "wilcox", only.pos = FALSE, logfc.threshold = 0.001, min.pct = 0.001)

# Save the output deg file
save(degs_of_subsets, file = "degs_of_cells_and_genes_subset.RData")

# Include gene IDs in the deg table
degs_of_subsets$gene_ID = rownames(degs_of_subsets)

# Find potential targets of RCO using the selection criteria defined in the function "target_chaser"
potential_markers = target_chaser(DEG_file = degs_of_subsets, store_dir = getwd(), store_outputs = FALSE)

# Lets save the output file
target_chaser(DEG_file = degs_of_subsets, store_dir = getwd(), store_outputs = TRUE, marker_file_name = "Direct_targets_of_RCO.csv")


