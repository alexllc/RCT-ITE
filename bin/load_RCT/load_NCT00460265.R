#' Script to load the SPECTRUM trial NCT00460265
#' Primary obj: overall survival (DTHDY, DTH)
#' Secondary obj: 
#' - Overall Objective Response Rate, (RDTLR, )
#' - Duration of Response, Time to Progression,
#' - Time to Response, 
#' - Progression Free Survival

source("./bin/load_lib.R")

file_path <- "./dat/PDS/HeadNe_Amgen_2007_265_NCT00460265/csv/"
corevar <- read.csv(paste0(file_path, "corevar.csv"))
eendpt <- read.csv(paste0(file_path, "aeendpt.csv"))
demo <- read.csv(paste0(file_path, "demo.csv"))

identical(corevar$STUDYID, eendpt$STUDYID) # check if the indexing is matched

# set to NA in R language for any type of missing or not otherwise specified entries
na_strings <- c("", "NA", "98", 98, "99", 99, "missing")
corevar <- corevar %>% replace_with_na_all(condition = ~.x %in% na_strings)
demo <- demo %>% replace_with_na_all(condition = ~.x %in% na_strings)

corevar_sel <- dplyr::select(corevar, all_of(c("SUBJID", "AGE", "SEX", "RACE", "PRHNTRTC", "TRTDUR", "B_ECOGCT", "DIAGTYCD", "TUMCAT", "DSTATUS", "PMAB", "HPVCD")))

# most the variable entries for EGFR related info is missing, no need to include
for (colid in 1:dim(demo)[2]) {
    print(colnames(demo)[colid])
    print(length(which(is.na(demo[,colid]))) / dim(demo)[1])
}

demo_sel <- dplyr::select(demo, all_of(c("SUBJID", "B_WEIGHT", "B_HEIGHT", "B_BSA", "PRSURG", "PRRADIO", "DIAGMONS", "METMONS", "DIAGSDCD", "HDIFFMCD", "HSSBTMCD", "DISTMET", "RECDIS", "CHILDPOT")))

X <- left_join(corevar_sel, demo_sel, by = c("SUBJID"))
X$SUBJID <- NULL
X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = as.data.frame(X), save_ddt = FALSE)

W <- as.numeric(corevar$ATRT == "panit. plus chemotherapy")

#
# Export trial data
#
NCT00460265 <- list(X_imp, W)

# outcome measurements
OS <- data.frame(T = eendpt$DTHDY / 30.4167, C = eendpt$DTH) # primary outcome
PFS <- data.frame(T = eendpt$PFSDYLR / 30.4167, C = eendpt$PFSLR) # secondary outcome
ORR <- data.frame(T = eendpt$RDYLR / 30.4167, C = eendpt$ROSLR) # secondary outcome

NCT00460265_outcomes <- c("OS", "PFS", "ORR")

for (outcome in NCT00460265_outcomes) {
    assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = get(outcome)[,2], X = X_imp)))
}
