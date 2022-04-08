source("./bin/load_lib.R")
source("./bin/impute_survival.R")

file_path <- "./dat/PDS/Breast_Allianc_2002_194_NCT00041119/csv/"

eval <- read.csv(paste0(file_path, "eval46_3_finala.csv"))

# set to NA in R language for any type of missing or not otherwise specified entries
missing_indicators <- c("", "NA", "98", 98, "99", 99)
for (missing in missing_indicators) {
    eval[eval == missing] <- NA
}

X <- dplyr::select(eval, all_of(c("RACE_ID", "stra1", "stra2", "OH002", "OH003", "OH004", "OH005", "OH006", "OH011", "OH012", "OH013", "OH014", "OH016", "OH027", "OH028", "OH032", "OH036", "OH037", "num_pos_nodes", "tsize", "agecat")))

X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

Y_list <- impute_survival(T = eval$dfsmos, C = eval$dfsstat, X = X_imp)

W <- as.numeric(eval$indrx == 2 | eval$indrx == 4) # Experimental arm is: 2=CA-6 or 4=T-6
NCT00041119_length <- list(X_imp, Y_list, W)

