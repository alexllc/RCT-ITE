
file_path <- "./dat/PDS/Breast_Allianc_2002_194_NCT00041119/csv/"

eval <- read.csv(paste0(file_path, "eval46_3_finala.csv"))
ae <- read.csv(paste0(file_path, "aeclean_finala.csv"))
# set to NA in R language for any type of missing or not otherwise specified entries
eval <- eval %>% replace_with_na_all(condition = ~.x %in% na_strings)
ae <- ae %>% replace_with_na_all(condition = ~.x %in% na_strings)

X <- dplyr::select(eval, all_of(c("RACE_ID", "stra1", "stra2", "OH002", "OH003", "OH004", "OH005", "OH006", "OH011", "OH012", "OH013", "OH014", "OH016", "OH027", "OH028", "OH032", "OH036", "OH037", "num_pos_nodes", "tsize", "agecat")))

X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

# Only compare CA vs Paclitaxel
W <- as.numeric(eval$indrx == 3 | eval$indrx == 4) # Experimental arm is: 3=T-4 or 4=T-6

NCT00041119_chemo <- list(X_imp, W)

toxicity <- ae %>% group_by(MASK_ID) %>% mutate(total_ae_grade = sum(GRADE_ID)) %>% select(all_of(c("MASK_ID", "total_ae_grade"))) %>% unique()
eval <- left_join(eval, toxicity, by = c("mask_id" = "MASK_ID"))


RFS <- data.frame(T = eval$dfsmos, C = eval$dfsstat)
dth <- eval$cod
dth[dth == 2] <- 1
OS <- data.frame(T = eval$survmos, C = dth)
AE <- eval$total_ae_grade

NCT00041119_chemo_outcomes <- c("RFS", "OS", "AE")

for (outcome in NCT00041119_chemo_outcomes) {
    if (outcome != "AE") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp)))
    } else {
        assign(paste0(outcome, "_Y_list"), list(get(outcome), NA, NA))
    }
}

save(NCT00041119_chemo, NCT00041119_chemo_outcomes, OS_Y_list, RFS_Y_list, AE_Y_list, file = "./bin/load_RCT/RCT_obj/NCT00041119_chemo.RData")