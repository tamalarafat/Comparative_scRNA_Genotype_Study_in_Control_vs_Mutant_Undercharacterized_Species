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

# Load the markers file
potential_markers = read.csv(file = "Potential_target_markers/Direct_targets_of_RCO.csv", row.names = 1)

for (j in c(1:nrow(potential_markers))){
  
  if (!dir.exists("Analysis_outputs")){
    dir.create("Analysis_outputs", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  if (!dir.exists(str_c("Analysis_outputs", "/", "Features_plot"))){
    dir.create(str_c("Analysis_outputs", "/", "Features_plot"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # Variable containing the directory path
  temp_dir = str_c("Analysis_outputs", "/", "Features_plot", "/")
  
  
  p <- FeaturePlot(integrated.data, 
                   features = rownames(potential_markers)[j], 
                   reduction = "umap", 
                   pt.size = 1, 
                   label = TRUE, 
                   label.size = 10,
                   split.by = "Genotypes", 
                   order = T, 
                   min.cutoff = 0.001, 
                   cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
    NoLegend() + 
    theme(
      line = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_blank())
  
  ggsave(plot = p, filename = str_c(temp_dir, "Feature_plot_", rownames(potential_markers)[j], ".png"), width = 8, height = 8, dpi = 300)
}
