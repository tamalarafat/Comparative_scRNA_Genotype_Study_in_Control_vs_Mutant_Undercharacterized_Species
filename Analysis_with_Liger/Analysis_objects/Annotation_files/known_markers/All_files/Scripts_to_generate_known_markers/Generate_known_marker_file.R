# Get the gene ids included in the data
col_genes = rownames(integrated.data)

# Known cell-type-specific markers
annotatedMarkers = read.xlsx("/netscratch/dep_tsiantis/grp_laurent/tamal/2021/Final_Analyses/Input_Files/Leaf-SC-marker-genes.xlsx")

####### The Epidermis markers ######
EpidermisMarkers = read.xlsx("/netscratch/dep_tsiantis/grp_laurent/tamal/2021/Final_Analyses/Input_Files/markers_epidermis.xlsx", colNames = T)
EpidermisMarkers <- EpidermisMarkers[c(10, 11, 12), c(2, 4, 5)]
colnames(EpidermisMarkers) <- c("ID", "Gene_Name", "Tissue")
rownames(EpidermisMarkers) <- NULL

rownames(EpidermisMarkers) <- EpidermisMarkers$ID

######### Epidermis markers till here

annotatedMarkers = annotatedMarkers[ , 1:3]
colnames(annotatedMarkers) <- c("ID", "Gene_Name", "Tissue")
rownames(annotatedMarkers) <- annotatedMarkers$ID
known_markers = rownames(annotatedMarkers) # Known cell-type-specific marker names

atids = intersect(col_genes, known_markers) # 34 known cell-type-specific markers are present in the OX WT data after filtering lowly expressed genes

KM_DF = annotatedMarkers[order(annotatedMarkers$Tissue, annotatedMarkers$Gene_Name), ] # Sorting the annotation file based on Tissue names

KM_names <- c("FIL", "KAN1", "KAN2", "YAB2", "YAB3", "AS2", "PHAVOLUTA", "PHB", "REVOLUTA", "CER2", "FAS2", "ATML1", "KCS10", "ANAC034", "AT1G11850", "BNQ2", 
              "LHCB8", "TMM", "LTP1", "CUT1", "PDF2", "GC1", "LHCA6", "LHCB7", "RBCS1B", "LHCB2.3", "EPFL9", "CA", "WOX1", "WOX3", "ATHB8", "PXY", "SMXL5", 
              "TMO6", "APL", "SWEET11", "TMO5", "VND7")

KM_DF$GeneName <- KM_names

# ,"BNQ2", "VND7", "FAS2","LHCB8", "APL"

HE = KM_DF[!KM_DF$GeneName %in% c("KAN2", "YAB3", "PHAVOLUTA", "CER2", "FAS2", "AT1G11850", "LHCB8", "GC1", "LHCB7", "WOX3", "APL", "VND7"), ]

KM_HE_DF <- HE[!HE$GeneName %in% c("ANAC034", "BNQ2", "FIL", "YAB2", "AS2", "PHB", "REVOLUTA", "WOX1", "APL"), ]
KM_HE_DF$Gene_Name <- KM_HE_DF$GeneName
KM_HE_DF$GeneName <- NULL
#### Lets create a dataframe according to Nora's markers order

CellTypeMarker <- rbind(KM_HE_DF, EpidermisMarkers)
CellTypeMarker <- CellTypeMarker[c(2, 7, 19, 20, 21, 8, 11, 12, 13, 15, 18, 16, 14), ]


# Creating a dataframe with the candidate genes only
RCO <- c("AT5G67651", "RCO", "NA")
STM <- c("AT1G62360", "STM", "NA")
KNAT1 <- c("AT4G08150", "KNAT1", "NA")
LMI1 <- c("AT5G03790", "LMI1", "NA")

Candidate_Genes <- data.frame(rbind(RCO, STM, KNAT1, LMI1))
colnames(Candidate_Genes) = c("AT_ID", "Gene_Name", "Tissue")
rownames(Candidate_Genes) <- Candidate_Genes$AT_ID

CellTypeMarker <- rbind(CellTypeMarker, Candidate_Genes)

# Load the new list of known cell-type markers
cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Thaliana_known_cell_type_markers.csv")
rownames(cell_type_markers) = cell_type_markers$AT_ID

# Get the ordered gene ids based on the cell-type
markers_order = CellTypeMarker$ID

# write_lines(markers_order, "known_cell_type_markers_ordered.txt")
known_markers = cell_type_markers
rownames(known_markers) = known_markers$AT_ID
known_markers = known_markers[markers_order, ]

present_genes = cell_type_markers[intersect(col_genes, rownames(cell_type_markers)), ]

present_genes = present_genes[c("AT1G77990", "AT5G41920", "AT1G46480", "AT5G65590", "AT1G12480", "AT5G25980", "AT3G01670", "AT3G01680", "AT1G08810", "AT2G46720", "AT1G79840", "AT5G41315", "AT3G27920"),]

curated_markers_file = rbind(known_markers[c(1:5), ], present_genes[c(4:6, 9:13), ], known_markers[c(6:8, 9:12), ], present_genes[c(1:3, 7, 8), ], known_markers[13, ], Candidate_Genes)

write.csv(curated_markers_file, "Curated_markers_thaliana.csv", row.names = FALSE)
