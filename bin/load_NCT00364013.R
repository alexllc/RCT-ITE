clin_df <- read.csv("./dat/NCT00364013_adls.csv")
biom <- read.csv("./dat/NCT00364013_biom.csv")

# Subjects with at least one biomarker entry is excluded from model building
tested_subj <- Reduce(union, lapply(biom[,c(2:dim(biom)[2])], function(X) which(is.na(X))))
model_clin_df <- clin_df[!clin_df$SUBJID %in% tested_subj,]

dim(model_clin_df)

model_clin_df$PFS_mo <- model_clin_df$PFSDYCR / 30.4167
model_clin_df$OS_mo <- model_clin_df$DTHDY / 30.4167

X <- select(model_clin_df, all_of(c("LIVERMET", "DIAGMONS", "AGE",  "SEX", "B_WEIGHT", "B_HEIGHT", "RACE",  "B_ECOG", "HISSUBTY", "B_METANM", "DIAGTYPE")))
W <- model_clin_df$ATRT

########################################################################
#### Step 2: address outcome censoring
########################################################################

# 1. Efron's tail correction
source("./audit/original_surv.r")
efron_tail <- exp(impute.survival(surv.time = model_clin_df$PFS_mo, censor = model_clin_df$PFSCR)) # no need to convert to log scale

# 2. Imputed tail correction (Khan & Shaw)
# sort matrix before imputing 
chron_os <- order(model_clin_df$PFS_mo)
imputeYn <- imputeYn(X = as.matrix(X[chron_os,]), Y = log(model_clin_df$PFS_mo[chron_os]), delta = model_clin_df$PFSCR[chron_os], method = "condMean")

# largest observation is not censored, no need to impute Yn

# use imputed Yn to perform Efron's tail correction again
efron_tail <- exp(impute.survival(surv.time = model_clin_df$PFS_mo, censor = model_clin_df$PFSCR))

# 3. Pseudo-observation
ps_mean <- pseudomean(time = model_clin_df$PFS_mo, event = model_clin_df$PFSCR, tmax = floor(max(model_clin_df$PFS_mo))) # negative values
