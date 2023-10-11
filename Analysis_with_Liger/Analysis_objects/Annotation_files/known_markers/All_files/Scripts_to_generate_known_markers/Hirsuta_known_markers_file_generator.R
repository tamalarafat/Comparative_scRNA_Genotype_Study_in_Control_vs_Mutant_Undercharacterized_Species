Known_markers_new = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/101122_Leaf_Cell_Type_markers_v12_with_ortho_ALL.csv")

Hirsuta_known_cell_type_markers = Known_markers_new[, c("Chirsutav12", "Short_name", "Tissue", "Athaliana")]
Hirsuta_known_cell_type_markers = na.omit(Hirsuta_known_cell_type_markers)
colnames(Hirsuta_known_cell_type_markers) = c("CH_ID", "Gene_Name", "Tissue", "AT_ID")

write.csv(Hirsuta_known_cell_type_markers, "Hirsuta_known_cell_type_markers.csv", row.names = FALSE)

# 22/02/2023
known_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Thaliana_known_cell_type_markers.csv")

cknown_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Hirsuta_known_cell_type_markers.csv")

target_markers = known_markers[c(83, 84, 85, 86), ]

# RCO, STM, KNAT1,  LMI1
target_genes = c("Chir06Ox-b35150.2", "Chir02Ox-b03100.2", "Chir08Ox-b10180.2", "Chir06Ox-b35160.2")
target_markers$CH_ID = target_genes

target_markers = target_markers[, c(4, 2, 3, 1)]
Hirsuta_known_cell_type_markers = rbind(cknown_markers, target_markers)

duplicated(df$column_name)


write.csv(Hirsuta_known_cell_type_markers, "/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Hirsuta_known_cell_type_markers_latest.csv", row.names = FALSE)


known_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Hirsuta_known_cell_type_markers_latest.csv")

Hirsuta_known_cell_type_markers = known_markers[!duplicated(known_markers$Gene_Name), ]
write.csv(Hirsuta_known_cell_type_markers, "/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Hirsuta_known_cell_type_markers_latest.csv", row.names = FALSE)

# 23/02/2023
cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Curated_markers_thaliana.csv")
cknown_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Hirsuta_known_cell_type_markers.csv")

intersect(cell_type_markers$AT_ID, cknown_markers$AT_ID)
a = cell_type_markers[cell_type_markers$AT_ID %in% intersect(cell_type_markers$AT_ID, cknown_markers$AT_ID), ]
b = merge(x = a, y = cknown_markers, by.x = "AT_ID", by.y = "AT_ID", all.y = FALSE)
b = b[!duplicated(b$AT_ID), c(1,2,3,4)]

write.csv(b, "/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Curated_markers_hirsuta.csv", row.names = FALSE)
