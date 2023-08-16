# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)


###
# Perform GO annotation of the genes in each GEP for different factorization
###

Thaliana_genes =  read.delim("/netscratch/dep_tsiantis/grp_laurent/tamal/2021/Final_Analyses/Input_Files/ATgenes.txt", header = FALSE, col.names = "genes")
Thaliana_genes <- as.character(Thaliana_genes$genes)

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1"

GO_annotation_of_GEPs(gene_clusters_file = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/Basis_objects/Basis.RData", 
                      ATgenes = Thaliana_genes, 
                      store_dir = storing_dir, 
                      store_folder = "GO_annotation_of_GEPs", 
                      use_ortho = TRUE, 
                      ortho_table = ortho_table, 
                      gene_column_in_use = "C.hirsutaOX", 
                      ortho_column_to_use = "A.thaliana.TAIR10")
