
source("./bin/impute_survival.R")

adsl <- read.csv("./dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/adsl_pds2019.csv")

# check if assigned trt is the same as actual trt
all(adsl$TRT==adsl$ATRT)
message(c("Which patient(s) had different assigned vs actual treatment: ", which(!(adsl$TRT==adsl$ATRT))))

# disgard assigned treatment
adsl$TRT <- NULL
adsl[adsl == ""] <- NA

# arm 1: Panitumumab + FOLFOX will be considered as W=1, arm 2: FOLFOX alone will be considered as W=0
adsl$ATRT <- as.numeric(adsl$ATRT == "Panitumumab + FOLFOX")
print(table(adsl$ATRT))

biomarker <- read.csv("./dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/biomark_pds2019.csv")
biom <- biomarker[,c(1,seq(3, dim(biomarker)[2], 2))]
biom[biom == ""] <- NA
for(col in colnames(biom)) {
    print(missing_prop <- sum(is.na(biom[[col]])) / dim(biom)[1])
}
# only BMMTR1 had sufficient non-missing data
clin_df <- left_join(adsl, as.data.frame(biom[,c(1:2)]), by = c("SUBJID"))

# convert categorical to numeric
num_df <- cat2num(clin_df)
num_clin_df <- num_df[[1]]
clin_df_ddt <- num_df[[2]]

#
# impute missing data if necessary
#
if (any(is.na(num_clin_df))) {
    pre_proc <- scale(num_clin_df, center = TRUE, scale = TRUE)
    imp_clin_df <- impute.knn(pre_proc, k=10)
    message(c("Seed used: ", imp_clin_df$rng.seed))
    imp_clin_df <- t(apply(imp_clin_df$data, 1, function(r) r *attr(imp_clin_df$data,'scaled:scale') + attr(imp_clin_df$data, 'scaled:center')))
}
message("Is there still missing data? ", any(is.na(imp_clin_df)))
imp_clin_df <- as.data.frame(imp_clin_df)
class(imp_clin_df$PFSCR) <- "integer" # after imputation returned floating point for 0s

# convert survival days to months
imp_clin_df$PFS_mo <- imp_clin_df$PFSDYCR / 30.4167
imp_clin_df$OS_mo <- imp_clin_df$DTHDY / 30.4167

# Format for causal inference

X <- select(imp_clin_df, all_of(c("LIVERMET", "DIAGMONS", "AGE",  "SEX", "B_WEIGHT", "B_HEIGHT", "RACE",  "B_ECOG", "HISSUBTY", "B_METANM", "DIAGTYPE", "BMMTR1")))
W <- adsl$ATRT[adsl$SUBJID %in% biom$SUBJID]
Y_list <- impute_survival(T = imp_clin_df$PFS_mo, C = imp_clin_df$PFSCR, X = X)

NCT00364013_KRASe2 <- list(X, Y_list, W)