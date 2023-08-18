# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF")

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Remove_CC_redo_clustering/Seurat_object_with_specified_UMAP/seurat_object_of_K_50_without_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

cluster.ids = levels(Idents(integrated.data))

for (i in c(1:length(cluster.ids))){
  heatmap_cells_coefficient(seuratObject = integrated.data, store_dir = getwd(), reduction_pattern = "inmfcc",store_folder = "Heatmap_of_clusters_inmfcc", cell_ids = WhichCells(integrated.data, idents = cluster.ids[i]), figureName = str_c("Cluster_", cluster.ids[i]))
}

# Lets plot the heatmaps with all the cells
map_factor_loadings(seuratObject = integrated.data, store_dir = getwd(), store_folder = "GEPs_distribution", factor_to_plot = "inmfcc")
