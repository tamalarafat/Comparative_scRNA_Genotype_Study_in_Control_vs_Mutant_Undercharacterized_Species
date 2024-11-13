# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_GEP20_22_and_CCGEPs/seurat_object_of_K_44_Without_GEP20_and_CC.RData")

# Storing directory
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_GEP20_22_and_CCGEPs/9_Trajectory_analysis"

############################ Change this chunk ############################
Idents(integrated.data) <- integrated.data$RNA_snn_res.0.3

paste("What is the active ident?", paste(levels(integrated.data@active.ident), collapse = ", "))

Idents(integrated.data) <- factor(Idents(integrated.data), levels = seq(0, length(levels(Idents(integrated.data))) - 1))
###########################################################################

############################# Trajectory Inference of the whole dataset
# While the slingshot vignette uses SingleCellExperiment, slingshot can also take a matrix of cell embeddings in reduced dimension as input. 
# We can optionally specify the cluster to start or end the trajectory based on biological knowledge.

# Here, I specified the STM cell cluster as the starting cluster for developmental progression.
# # This way we learn the cluster relationships in a semi-supervised manner by incorporating prior knowledge of established gene marker (STM)
Species_TI <- slingshot(Embeddings(integrated.data, "inmfcc"), clusterLabels = Idents(integrated.data), start.clus = 13, stretch = 0, maxit = 10)

save(Species_TI, file = "WT_CH_TI.RData")
