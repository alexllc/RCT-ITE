#' Load dataset for NCT00052910
#' 
#' Pirmary outcome: Overall survival
#' Secondary outcome: Progression free survival
#' 

file_path <- "./dat/PDS/Multiple_Allianc_2002_213_NCT00052910/"

eff <- read.csv(paste0(file_path, "NCT00052910_D1_(EFFICACY).csv"))
trial_na_strings <- c("", "NA","<NA>", ".", "-2", -2, "-1", -1)
eff <- eff %>% replace_with_na_all(condition = ~.x %in% trial_na_strings)

# fix specific variable missings
eff$PD_location[eff$PD_location == 0] <- NA
eff$ETHNIC_ID[eff$ETHNIC_ID == 9] <- NA

X <- dplyr::select(eff, -all_of(c("MASK_ID", "TREAT_ASSIGNED", "death_final", "os_year", "DFS_status", "DFS_year")))
X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = as.data.frame(X), save_ddt = FALSE)

# treatment assignments
W <- as.numeric(eff$TREAT_ASSIGNED == 1)
NCT00052910<- list(X_imp, W)

# outcome measurements

OS <- data.frame(T = eff$os_year, C = eff$death_final) # primary outcome
PFS <- data.frame(T = eff$DFS_year, C = eff$DFS_status) # secondary 
NCT00052910_outcomes <- c("OS", "PFS")

for (outcome in NCT00052910_outcomes) {
    assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp)))
}

save(NCT00052910, NCT00052910_outcomes, OS_Y_list, PFS_Y_list, file = "./bin/load_RCT/RCT_obj/NCT00052910.RData")