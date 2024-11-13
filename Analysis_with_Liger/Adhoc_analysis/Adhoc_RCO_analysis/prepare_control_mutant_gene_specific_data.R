# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_rco_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

integrated.data$cell_ID = rownames(integrated.data@meta.data)

# Check the RCO expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", ] != 0])

# Subset the seurat object
rco.cells.subset = subset(integrated.data, subset = cell_ID %in% Cells_with_expression_detection)

save(rco.cells.subset, file = "rco_cells_subset_wt_mutant.RData")  

# Genes to test - combine - GEP32, cluster 5 (liger), cluster 7 (seurat) marker genes

GEP32 = read.delim("GEP_32_genes.txt", header = FALSE)[, 1]
SM_cl5_liger = read.delim("Cluster_5_specific_markers_CH_ID.txt", header = FALSE)[, 1]
SM_cl7_seurat = read.csv("seurat_rco_cluster_7.csv")[, 8]

intersect(GEP32, SM_cl7_seurat) # 22 genes
intersect(GEP32, SM_cl5_liger) # 23 genes

combined_rco_markers = unique(c(GEP32, SM_cl5_liger, SM_cl7_seurat)) # 124 genes

writeLines(combined_rco_markers, "all_rco_markers.txt")

