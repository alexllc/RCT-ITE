source("./bin/load_lib.R")
source("./bin/impute_survival.R")

file_path <- "./dat/PDS/LungSm_Allianc_1998_261_NCT00003299/csv/"

demo <- read.csv(paste0(file_path, "c9732_demographic.csv"))
demo <- filter(demo, ELIGIBLE == 2)

W <- demo$TRT_ARM - 1

X <- dplyr::select(demo, all_of(c("GENDER", "AGE", "RACE", "PS", "WGT_LOSSPCT", "NUM_META", "SCLC_STAGE", "CHEMO_CYCLE", "RADIOTHERAPY", "OFFTRT_RX")))

X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

Y_list <- impute_survival(T = demo$PFS_TIME, C = demo$PFS_STATUS, X = X_imp)

NCT00003299 <- list(X_imp, Y_list, W)