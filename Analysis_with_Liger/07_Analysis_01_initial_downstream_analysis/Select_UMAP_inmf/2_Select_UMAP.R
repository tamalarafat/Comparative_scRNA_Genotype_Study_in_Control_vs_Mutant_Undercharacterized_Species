# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_object/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

############################ Set the resolution parameter ############################

Idents(integrated.data) <- factor(Idents(integrated.data), levels = seq(0, length(levels(Idents(integrated.data))) - 1))

paste("What is the active ident?", paste(levels(integrated.data@active.ident), collapse = ", "))

###########################################################################

integrated.data <- RunUMAP(integrated.data, reduction = "inmfcc", dims = 1:ncol(integrated.data@reductions$inmf@cell.embeddings), n.components = 2, n.neighbors = 40, seed.use = 10)

save(integrated.data, file = "seurat_object_of_K_50.RData")
