source("./bin/impute_survival.R")


corevar <- read.csv("./dat/PDS/Colorec_Amgen_2005_262_NCT00115765/csv/corevar.csv")
eendpt <- read.csv("./dat/PDS/Colorec_Amgen_2005_262_NCT00115765/csv/a_eendpt.csv")
identical(corevar$STUDYID, eendpt$STUDYID) # check if the indexing is matched
corevar[corevar == "" | corevar == "NA"] <- NA

X <- dplyr::select(corevar, all_of(c("AGE", "SEXCD", "RACCATCD", "ACHEMOCD", "B_LESNM", "B_LESNMT", "B_LESNMN", "B_METANM", "FSTOXDS", "DURSURG", "PRADJICD", "B_LDHYN", "B_ECOGCD", "DIAGTYCD", "PDADEN", "KRASCD", "TRTDUR")))
X <- missing_too_much(X)


W <- as.numeric(corevar$ATRT == "panit. plus bevacizumab with chemotherapy")

X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

Y_list <- impute_survival(T = eendpt$PFSDYCR / 30.4167, C = eendpt$PFSCR, X = X_imp)

NCT00115765 <- list(X_imp, Y_list, W)