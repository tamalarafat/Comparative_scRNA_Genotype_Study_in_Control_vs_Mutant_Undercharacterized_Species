# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Remove_CC_redo_clustering/seurat_object_of_K_50_without_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

############################ Set the resolution parameter ############################

Idents(integrated.data) <- factor(Idents(integrated.data), levels = seq(0, length(levels(Idents(integrated.data))) - 1))

paste("What is the active ident?", paste(levels(integrated.data@active.ident), collapse = ", "))

###########################################################################

# Storing directory
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Remove_CC_redo_clustering"

# Lets get the metadata file
md = integrated.data@meta.data

# UMAP - 2 components reduction
n_neighbors = c(20, 30, 40, 50)
n_seed = c(1:50)

for (i in c(1:length(n_neighbors))){
  for (j in c(1:length(n_seed))){
    integrated.data <- RunUMAP(integrated.data, reduction = "inmfcc", dims = 1:ncol(integrated.data@reductions$inmfcc@cell.embeddings), n.components = 2, n.neighbors = n_neighbors[i], seed.use = n_seed[j])
    
    RDimension_plot(seuratObject = integrated.data, store_dir = res_dir, store_folder = "Seed_and_NN",dimension_reduction_name = "umap", figure_name_pref = str_c("_nn", n_neighbors[i], "_seed", n_seed[j]))
    }
}


