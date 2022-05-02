source("./bin/impute_survival.R")

adsl <- read.csv("./dat/PDS/Colorec_Amgen_2004_310_NCT00113763/csv/adsl_pds2019.csv")

biomark <- read.csv("./dat/PDS/Colorec_Amgen_2004_310_NCT00113763/csv/biomark_pds2019.csv")

biomarknm <- unlist(biomark[1,seq(2, dim(biomark)[2], 2)])
biomarknm <- gsub(" ", "_", biomarknm)

biom <- biomark[,c(1,seq(3, dim(biomark)[2], 2))]
colnames(biom) <- c("SUBJID", biomarknm)
biom[biom == ""] <- NA

clin_df = biom
for(col in colnames(clin_df)) {
    print(missing_prop <- sum(is.na(clin_df[[col]])) / dim(clin_df)[1])
}

identical(adsl$SUBJID, biom$SUBJID)

adsl <- left_join(adsl, biom, by = c("SUBJID"))

# check if assigned trt is the same as actual trt
all(adsl$TRT==adsl$ATRT)
message(c("Which patient(s) had different assigned vs actual treatment: ", which(!(adsl$TRT==adsl$ATRT))))
adsl[adsl == ""] <- NA

X <- dplyr::select(adsl, -all_of(c("SUBJID", "TRT", "ATRT", "DTHDYX", "DTHX", "PFSDYCR", "PFSCR")))

X <- missing_too_much(X)

X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

Y_list <- impute_survival(T = adsl$PFSDYCR / 30.4167, C = adsl$PFSCR, X = X_imp)

W <- as.numeric(adsl$ATRT == "panit. plus best supportive care")
print(table(W))

NCT00113763 <- list(X_imp, Y_list, W)