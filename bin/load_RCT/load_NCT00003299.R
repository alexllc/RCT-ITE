
file_path <- "./dat/PDS/LungSm_Allianc_1998_261_NCT00003299/csv/"

demo <- read.csv(paste0(file_path, "c9732_demographic.csv"))
demo <- demo %>% filter(ELIGIBLE == 2) %>% replace_with_na_all(condition = ~.x %in% na_strings)

ae <- read.csv(paste0(file_path, "c9732_ae.csv"))
ae <- ae %>% replace_with_na_all(condition = ~.x %in% na_strings)

W <- demo$TRT_ARM - 1

X <- dplyr::select(demo, all_of(c("GENDER", "AGE", "RACE", "PS", "WGT_LOSSPCT", "NUM_META", "SCLC_STAGE", "CHEMO_CYCLE", "RADIOTHERAPY")))

X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

NCT00003299 <- list(X_imp, W)

toxicity <- ae %>% group_by(PHATOM_ID) %>% mutate(total_ae_grade = sum(AE_GRADE)) %>% select(all_of(c("PHATOM_ID", "total_ae_grade"))) %>% unique()
demo <- left_join(demo, toxicity, by = c("PHATOM_ID"))

imp_outcome <- dplyr::select(demo, all_of(c("PFS_TIME", "PFS_STATUS", "PD_TIME", "PD", "OS_TIME", "STATUS", "total_ae_grade")))
imp_outcome$PD <- imp_outcome$PD - 1
imp_outcome$STATUS <- imp_outcome$STATUS - 1
imp_outcome$STATUS[imp_outcome$STATUS == 2] <- 0
imp_outcome <- impute_df_missing(clin_df = as.data.frame(imp_outcome), save_ddt = FALSE)

PFS <- data.frame(T = imp_outcome$PFS_TIME, C = ceiling(imp_outcome$PFS_STATUS))
PD <- data.frame(T = imp_outcome$PD_TIME, C = ceiling(imp_outcome$PD))
OS <- data.frame(T = imp_outcome$OS_TIME, C = ceiling(imp_outcome$STATUS))
AE <- imp_outcome$total_ae_grade

NCT00003299_outcomes <- c("PFS", "PD", "OS", "AE")

for (outcome in NCT00003299_outcomes) {
    if (outcome != "AE") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp)))
    } else {
        assign(paste0(outcome, "_Y_list"), list(get(outcome), NA, NA))
    }
}