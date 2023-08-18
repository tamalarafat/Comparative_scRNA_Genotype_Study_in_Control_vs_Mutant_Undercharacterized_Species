# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_leaf_protoplast_v12_final_table_DEGs_2reps_final_August_2022.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# WT C. hirsuta
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
OX_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_OX_RNA_1ST_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
OX_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_2nd_ALL_2_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 3rd Experiment (leaf 5 and 6)
OX_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_3rd_ALL_3000_Newest/filtered_feature_bc_matrix/")

# Load data - WT OX 7th Experiment (leaf 6 and 7)
OX_data_7E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_7th_ALL_2_Newest/filtered_feature_bc_matrix/")

###
# rco C. hirsuta
###

# WT OX 3rd Experiment - leaf 5 and 6
rco_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_rco_RNA_3rd_ALL_Newest/filtered_feature_bc_matrix/")

# WT OX 3rd Experiment - leaf 5 and 6
rco_data_6E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_rco_RNA_6th_Newest/filtered_feature_bc_matrix/")

# WT OX 7th Experiment - leaf 6 and 7
rco_data_7E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_rco_RNA_7th_ALL_2_Newest/filtered_feature_bc_matrix/")

# All gene IDs - Cardamine hirsuta
genes_hirsuta = rownames(OX_data_1E)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(genes_hirsuta, PP_genes)

# Subset the data tables - without protoplasting-induced genes
OX_data_1E <- OX_data_1E[genes_to_keep, ]
OX_data_2E <- OX_data_2E[genes_to_keep, ]
OX_data_3E <- OX_data_3E[genes_to_keep, ]
OX_data_7E <- OX_data_7E[genes_to_keep, ]

rco_data_3E <- rco_data_3E[genes_to_keep, ]
rco_data_6E <- rco_data_6E[genes_to_keep, ]
rco_data_7E <- rco_data_7E[genes_to_keep, ]

# Create seurat object and perform initial filtering - 
# 1. Remove genes if their expression was not detected in at least one cell out of every 500 ("min.cells"),
# 2. Remove cells if at least 200 genes were not detected to be expressed (min.features = 200),
# 3. Remove cells with a total count of more than 110000 (nCount_RNA > 110000).
# 4. Remove cells if 5% or more of the total count of a cell belongs to the mitochondiral genes.
# 5. Remove cells if 10% or more of the total count of a cell belongs to the chloroplast genes.

###
# OX - 1 E
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_1E <- CreateSeuratObject(counts = OX_data_1E, project = "OX_1E", min.cells = 13, min.features = 200)

# Add metadata information to the seurat object
OX_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_1E <- subset(OX_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_1E[["percent.mt"]] <- PercentageFeatureSet(OX_1E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_1E[["percent.pt"]] <- PercentageFeatureSet(OX_1E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_1E <- subset(OX_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_1E <- NormalizeData(OX_1E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
OX_1E <- FindVariableFeatures(OX_1E, selection.method = "vst", nfeatures = 2000)


###
# OX - 2 E
###

# First replicate - OX 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
OX_2E <- CreateSeuratObject(counts = OX_data_2E, project = "OX_2E", min.cells = 21, min.features = 200)

# Add metadata information to the seurat object
OX_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_2E <- subset(OX_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_2E[["percent.mt"]] <- PercentageFeatureSet(OX_2E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_2E[["percent.pt"]] <- PercentageFeatureSet(OX_2E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_2E <- subset(OX_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_2E <- NormalizeData(OX_2E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
OX_2E <- FindVariableFeatures(OX_2E, selection.method = "vst", nfeatures = 2000)


###
# OX - 3 E
###

# First replicate - OX 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
OX_3E <- CreateSeuratObject(counts = OX_data_3E, project = "OX_3E", min.cells = 8, min.features = 200)

# Add metadata information to the seurat object
OX_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_3E <- subset(OX_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_3E[["percent.mt"]] <- PercentageFeatureSet(OX_3E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_3E[["percent.pt"]] <- PercentageFeatureSet(OX_3E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_3E <- subset(OX_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_3E <- NormalizeData(OX_3E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
OX_3E <- FindVariableFeatures(OX_3E, selection.method = "vst", nfeatures = 2000)


###
# OX - 7 E
###

# First replicate - OX 7E - total cells 9090; filter out genes that are not detected in at least 18 cells
OX_7E <- CreateSeuratObject(counts = OX_data_7E, project = "OX_7E", min.cells = 18, min.features = 200)

# Add metadata information to the seurat object
OX_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-7", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_7E <- subset(OX_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_7E[["percent.mt"]] <- PercentageFeatureSet(OX_7E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_7E[["percent.pt"]] <- PercentageFeatureSet(OX_7E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_7E <- subset(OX_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_7E <- NormalizeData(OX_7E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
OX_7E <- FindVariableFeatures(OX_7E, selection.method = "vst", nfeatures = 2000)


###
# rco - 3 E
###

# First replicate - rco 3E - total cells 4000; filter out genes that are not detected in at least 8 cells
rco_3E <- CreateSeuratObject(counts = rco_data_3E, project = "rco_3E", min.cells = 8, min.features = 200)

# Add metadata information to the seurat object
rco_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "rco-OX-3", "rco", "Leaf")

# Remove cells with a total count more than 110000
rco_3E <- subset(rco_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
rco_3E[["percent.mt"]] <- PercentageFeatureSet(rco_3E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
rco_3E[["percent.pt"]] <- PercentageFeatureSet(rco_3E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
rco_3E <- subset(rco_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
rco_3E <- NormalizeData(rco_3E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
rco_3E <- FindVariableFeatures(rco_3E, selection.method = "vst", nfeatures = 2000)


###
# rco - 6 E
###

# First replicate - rco 6E - total cells 9200; filter out genes that are not detected in at least 18 cells
rco_6E <- CreateSeuratObject(counts = rco_data_6E, project = "rco_6E", min.cells = 18, min.features = 200)

# Add metadata information to the seurat object
rco_6E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "rco-OX-6", "rco", "Leaf")

# Remove cells with a total count more than 110000
rco_6E <- subset(rco_6E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
rco_6E[["percent.mt"]] <- PercentageFeatureSet(rco_6E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
rco_6E[["percent.pt"]] <- PercentageFeatureSet(rco_6E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
rco_6E <- subset(rco_6E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
rco_6E <- NormalizeData(rco_6E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
rco_6E <- FindVariableFeatures(rco_6E, selection.method = "vst", nfeatures = 2000)


###
# rco - 7 E
###

# First replicate - rco 7E - total cells 5860; filter out genes that are not detected in at least 12 cells
rco_7E <- CreateSeuratObject(counts = rco_data_7E, project = "rco_7E", min.cells = 12, min.features = 200)

# Add metadata information to the seurat object
rco_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "rco-OX-7", "rco", "Leaf")

# Remove cells with a total count more than 110000
rco_7E <- subset(rco_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
rco_7E[["percent.mt"]] <- PercentageFeatureSet(rco_7E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
rco_7E[["percent.pt"]] <- PercentageFeatureSet(rco_7E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
rco_7E <- subset(rco_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
rco_7E <- NormalizeData(rco_7E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
rco_7E <- FindVariableFeatures(rco_7E, selection.method = "vst", nfeatures = 2000)

###
# Identifying features, anchors, and performing integration
###

# Integration of the replicates - find anchors
anchFeatures <- SelectIntegrationFeatures(object.list = list(OX_1E, OX_2E, OX_3E, OX_7E, rco_3E, rco_6E, rco_7E))

ingAnchors <- FindIntegrationAnchors(object.list = list(OX_1E, OX_2E, OX_3E, OX_7E, rco_3E, rco_6E, rco_7E), dims = 1:50, anchor.features = anchFeatures)

# To keep records of all the genes in the integrated assay, create a feature set with all of the genes from different replicates.
features_integrated <- unique(c(rownames(OX_1E), rownames(OX_2E), rownames(OX_3E), rownames(OX_7E), rownames(rco_3E), rownames(rco_6E), rownames(rco_7E)))

# Integrate the replicates
integrated.data <- IntegrateData(anchorset = ingAnchors, dims = 1:50, verbose = T, features.to.integrate = features_integrated)

# Setting the default assay to "integrated"
DefaultAssay(integrated.data) <- "integrated"

# Gene level scaling - standardization
integrated.data <- ScaleData(integrated.data, verbose = FALSE)

# Run PCA
integrated.data <- RunPCA(integrated.data, npcs = 50, verbose = FALSE)

# Run UMAP and tSNE
integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = 1:50, n.components = 3)

integrated.data <- RunTSNE(integrated.data, reduction = "pca", dims = 1:50, dim.embed = 2)

# Find neighbours and clusters
integrated.data <- FindNeighbors(integrated.data, reduction = "pca", dims = 1:50)

for (i in seq(0.1, 1.2, 0.1)) {
  integrated.data <- FindClusters(integrated.data, resolution = i, n.start = 50, n.iter = 50)
}

save(integrated.data, file = "integrated_ox_wt_rco_seurat.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_ox_wt_rco_seurat.txt")