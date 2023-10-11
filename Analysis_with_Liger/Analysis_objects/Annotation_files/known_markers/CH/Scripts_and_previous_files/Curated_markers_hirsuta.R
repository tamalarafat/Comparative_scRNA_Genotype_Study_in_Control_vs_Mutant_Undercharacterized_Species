# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Load the ortho table to add some of target genes' hirsuta IDs
ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/comparative_study_of_two_species/Liger_Analyses/Liger_Analysis_strategy_2/On_coefficient/Seurat_objects/seurat_object_of_K_50.RData")

# Cell type markers file for A. thaliana
cell_type_markers_AT = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Curated_markers_thaliana.csv")

# Assign row names to Arabidopsis IDs
rownames(cell_type_markers_AT) = cell_type_markers_AT$AT_ID

# Known cell type markers
cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Curated_markers_hirsuta.csv")

# Change the column names that is the same as thaliana known cell type marker file 
colnames(cell_type_markers) = c("AT_ID", "Gene_Name", "Tissue", "CH_ID")

# Get the subset of traget gene IDs
targets = tail(cell_type_markers_AT, 4)

# Get the hirsuta ids for the target genes
targets$CH_ID = c("Chir06Ox-b35150.2", "Chir02Ox-b03100.2", "Chir08Ox-b10180.2", "Chir06Ox-b35160.2")

# Add the rows, the target genes to hirsuta known cell-type markers
cell_type_markers = rbind.data.frame(cell_type_markers, targets)

rownames(cell_type_markers) = cell_type_markers$AT_ID

# A final check by looking at the known cell type markers 
known_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Hirsuta_known_cell_type_markers_latest.csv")

known_markers_subset = known_markers[known_markers$AT_ID %in% cell_type_markers$AT_ID, ]

# Yes, the gene ids, thaliana and hirsuta, are the correct ones

# Lets order the gene IDs of the dataframe by tissue
cell_type_markers = cell_type_markers[rownames(cell_type_markers_AT), ]

rownames(cell_type_markers) <- NULL

write.csv(cell_type_markers, "Curated_markers_hirsuta.csv", row.names = FALSE)
