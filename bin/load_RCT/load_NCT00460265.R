source("./bin/load_lib.R")
source("./bin/impute_survival.R")

file_path <- "./dat/PDS/HeadNe_Amgen_2007_265_NCT00460265/csv/"
corevar <- read.csv(paste0(file_path, "corevar.csv"))
eendpt <- read.csv(paste0(file_path, "aeendpt.csv"))
demo <- read.csv(paste0(file_path, "demo.csv"))

# OTHCANTX - other cancer tx for this cancer type
# RADIOTX - radiotherapy, mainly for the prior radiotx column


identical(corevar$STUDYID, eendpt$STUDYID) # check if the indexing is matched

# set to NA in R language for any type of missing or not otherwise specified entries
missing_indicators <- c("", "NA", "98", 98, "99", 99)
corevar[corevar %in% missing_indicators] <- NA
demo[corevar %in% missing_indicators] <- NA

corevar_sel <- dplyr::select(corevar, all_of(c("SUBJID", "AGE", "SEX", "RACE", "PRHNTRTC", "B_ECOGCT", "DIAGTYCD", "PMAB", "HPVCD")))

# most the variable entries for EGFR related info is missing, no need to include
for (colid in 54:78) {
    print(length(which(is.na(demo[,colid]))) / dim(demo)[1])
}

demo_sel <- dplyr::select(demo, all_of(c("SUBJID", "B_WEIGHT", "B_HEIGHT", "B_BSA", "PRSURG", "PRRADIO", "DIAGSDCD")))

X <- left_join(corevar_sel, demo_sel, by = c("SUBJID"))
X$SUBJID <- NULL

X <- missing_too_much(X)


W <- as.numeric(corevar$ATRT == "panit. plus chemotherapy")

X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

Y_list <- impute_survival(T = eendpt$PFSDYLR / 30.4167, C = eendpt$PFSLR, X = X_imp)

NCT00460265 <- list(X_imp, Y_list, W)
