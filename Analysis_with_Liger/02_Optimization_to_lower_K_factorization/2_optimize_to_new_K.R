# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Optimize to a new K from previously factorized K (Higher than the new K) 
###

storing_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1"

optimize_to_new_k(optimized_liger_object = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/comparative_study_of_wt_and_mutant_genotype_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/1_Integration_script/integrated_ox_wt_rco_liger.RData", 
                  choice_of_new_K = c(10:49),
                  store_dir = storing_dir, 
                  store_folder = "Liger_objects")
