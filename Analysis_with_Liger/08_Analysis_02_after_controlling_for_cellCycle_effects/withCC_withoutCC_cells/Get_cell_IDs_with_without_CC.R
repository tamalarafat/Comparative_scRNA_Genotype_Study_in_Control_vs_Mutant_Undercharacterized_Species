# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_3/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# With cell-cycle GEPs

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_3/Analysis_objects/Seurat_object_with_CC/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

md_with = integrated.data@meta.data
md_with$cell_ID = rownames(md_with)

md_with = md_with[, c("cell_ID", "RNA_snn_res.0.3")]
colnames(md_with)[2] <- "withCC"

save(md_with, file = "withCC.RData")

# Let's save the metadata information too
save(md_with, file = "metadata_withCC.RData")

# Without cell-cycle GEPs

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_3/Analysis_objects/Seurat_object_with_specified_UMAP/seurat_object_of_K_50_without_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

md_withoutCC = integrated.data@meta.data
md_withoutCC$cell_ID = rownames(md_withoutCC)

md_withoutCC = md_withoutCC[, c("cell_ID", "RNA_snn_res.0.3")]
colnames(md_withoutCC)[2] <- "withoutCC"

save(md_withoutCC, file = "withoutCC.RData")

# Let's save the metadata information too
save(md_withoutCC, file = "metadata_withoutCC.RData")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_wt_and_mutant_genotype_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_rco_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

integrated.data$reps <- integrated.data$orig.ident
integrated.data$reps <- factor(integrated.data$reps, levels = c("OX_1E", "OX_2E", "OX_3E", "OX_7E", "rco_3E", "rco_6E", "rco_7E"), labels = c("O1", "O2", "O3", "O7", "M3", "M6", "M7"))

md_with_seurat = integrated.data@meta.data
md_with_seurat$cell_ID = rownames(md_with_seurat)

# Let's save the metadata information too
save(md_with_seurat, file = "metadata_withCC_seurat.RData")

md_with_seurat = md_with_seurat[, c("cell_ID", "integrated_snn_res.0.4", "reps")]

# Add string to compare
md_with_seurat$cell_ID = str_c(md_with_seurat$reps, "_", substr(md_with_seurat$cell_ID, start = 1, stop = 18))

# get the columns of interest
md_with_seurat = md_with_seurat[, c(1, 2)]

colnames(md_with_seurat)[2] <- "withCC"

save(md_with_seurat, file = "withCC_seurat.RData")



