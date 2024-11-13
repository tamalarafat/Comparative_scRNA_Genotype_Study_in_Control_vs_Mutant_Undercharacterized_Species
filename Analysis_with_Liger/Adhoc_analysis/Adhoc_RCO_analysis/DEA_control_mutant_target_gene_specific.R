# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the RCO cells marker - chapter 1
rco_markers = read.delim("all_rco_markers.txt", header = FALSE)[, 1]

# load the seurat object
integrated.rco <- loadRData("rco_cells_subset_wt_mutant.RData")

integrated.rco$Genotype <- factor(integrated.rco$Genotype)

Idents(integrated.rco) <- "Genotype"

DefaultAssay(integrated.rco) <- "RNA"

rco_WT_markers = FindMarkers(integrated.rco, ident.1 = "WT", features = rco_markers, test.use = "wilcox", only.pos = FALSE, logfc.threshold = 0.001, min.pct = 0.001)
save(rco_WT_markers, file = "rco_WT_markers.RData")

rco_WT_markers$gene_ID = rownames(rco_WT_markers)

rco_WT_markers = rco_WT_markers[(rco_WT_markers$pct.1 - rco_WT_markers$pct.2) >= rco_WT_markers$pct.2, ]
save(rco_WT_markers, file = "rco_WT_markers_final_set.RData")

# Load the gene description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_3/Home_PC/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

cluster_markers_des = ATH[rownames(ATH) %in% rownames(rco_WT_markers), ]
cluster_markers_des = cluster_markers_des[, c(1, 2, 7)]

markers = merge(rco_WT_markers, cluster_markers_des, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[,c(1, 7, 8)]
rownames(markers) = markers$gene_ID
colnames(markers) = c("CH_ID", "AT_ID", "Gene_name")

save(markers, file = "rco_WT_markers_final_set_details.RData")

for (j in c(1:nrow(markers))){
  p <- FeaturePlot(integrated.rco, features = rownames(markers)[j], reduction = "umap", pt.size = 1, split.by = "Genotype", order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
    NoLegend() + 
    theme(
      line = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_blank())
  
  ggsave(plot = p, filename = if (!is.na(markers[j, "Gene_name"])){str_c("Feature_plot_", markers[j, "Gene_name"], ".png")} else {str_c("Feature_plot_", rownames(markers)[j], ".png")}, width = 6, height = 6, dpi = 400)
}
