projects_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# projects_dir = "~/Documents/Projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_1 <- list.files(paste0(projects_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE)
sapply(scripts_list_1, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(projects_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE)
sapply(scripts_list_2, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the file containing grouping of the genes from the basis matrix for each factorization
load(paste0(projects_dir, "comparative_study_of_wt_and_mutant_genotype_Cardamine/Analysis_with_Liger/Analysis_objects/Gene_clusters/Basis.RData"))

# Load the orthologues table
ortho_table = read.csv(paste0(projects_dir, "comparative_study_of_wt_and_mutant_genotype_Cardamine/Analysis_with_Liger/Analysis_objects/Annotation_files/Orthologues_table/Orthos_table.csv"))

# The gene IDs of interest has to be stored in the gene_ID column of the orthologues table
ortho_table$gene_ID = str_replace(string = ortho_table$C.hirsutaOX, pattern = "_", replacement = "-")

storing_dir = paste0(projects_dir, "comparative_study_of_wt_and_mutant_genotype_Cardamine/Analysis_with_Liger/Analysis_output")

Basis_to_GEP_generator(input_data = Basis, 
                       Factor_ID = "50", 
                       store_GEPs = TRUE, 
                       return_GEPs = FALSE, 
                       generate_ortho = FALSE, 
                       store_dir = storing_dir)

Basis_to_GEP_generator(input_data = Basis, 
                       Factor_ID = "50", 
                       store_GEPs = TRUE, 
                       return_GEPs = FALSE, 
                       generate_ortho = TRUE, 
                       ortho_file = ortho_table, 
                       conversion_ID = "A.thaliana.TAIR10", 
                       store_dir = storing_dir)

