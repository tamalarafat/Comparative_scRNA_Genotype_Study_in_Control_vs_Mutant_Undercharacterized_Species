Known_markers_new = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/101122_Leaf_Cell_Type_markers_v12_with_ortho_ALL.csv")

Thaliana_known_cell_type_markers = Known_markers_new[, c("Athaliana", "Short_name", "Tissue")]
Thaliana_known_cell_type_markers$Athaliana[44] = "AT1G12090"
colnames(Thaliana_known_cell_type_markers) = c("AT_ID", "Gene_Name", "Tissue")
Thaliana_known_cell_type_markers = Thaliana_known_cell_type_markers[-c(40, 57, 84, 86, 87, 88), ]

# Add these traget genes to the markers file
RCO <- c("AT5G67651", "RCO", "NA")
STM <- c("AT1G62360", "STM", "NA")
KNAT1 <- c("AT4G08150", "KNAT1", "NA")
LMI1 <- c("AT5G03790", "LMI1", "NA")

Candidate_Genes <- data.frame(rbind(RCO, STM, KNAT1, LMI1), row.names = NULL)
colnames(Candidate_Genes) = c("AT_ID", "Gene_Name", "Tissue")

write.csv(Candidate_Genes, "Thaliana_candidate_genes.csv", row.names = FALSE)

Thaliana_known_cell_type_markers <- rbind(Thaliana_known_cell_type_markers, Candidate_Genes)
write.csv(Thaliana_known_cell_type_markers, "Thaliana_known_cell_type_markers.csv", row.names = FALSE)
