
########################################################################
#### Step 1: Prepare data
########################################################################
clindat <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Multiple_Allianc_2002_213_NCT00052910/NCT00052910_D1_(EFFICACY).csv")
tx <- clindat$TREAT_ASSIGNED - 1
outcome_col <- c("death_final", "os_year", "DFS_status", "DFS_year")
outcome <- dplyr::select(clindat, all_of(outcome_col))
covar <- as.data.frame(dplyr::select(clindat, -c("MASK_ID", "TREAT_ASSIGNED", outcome_col)))

# prepare imputed covariate matrix (imported from Python)
covar <- read.csv("./audit/imputed_covar.csv")
covar$os_year <- NULL
covar$TREAT_ASSIGNED <- NULL

########################################################################
#### Step 2: address outcome censoring
########################################################################

# 1. Efron's tail correction
source("./audit/original_surv.r")
efron_tail <- exp(impute.survival(surv.time = clindat$os_year, censor = clindat$death_final)) # no need to convert to log scale

# 2. Imputed tail correction (Khan & Shaw)
# sort matrix before imputing 
chron_os <- order(clindat$os_year)
imputeYn <- imputeYn(X = as.matrix(covar[chron_os,]), Y = log(clindat$os_year[chron_os]), delta = clindat$death_final[chron_os], method = "condMean")
imputedYn <- exp(imputeYn$Yn)

# use imputed Yn to perform Efron's tail correction again
replace_Yn_os <- clindat$os_year
replace_Yn_os[order(clindat$os_year, decreasing = T)[1]] <- imputedYn
efron_tail_Yn <- exp(impute.survival(surv.time = replace_Yn_os, censor = clindat$death_final))

# 3. Pseudo-observation
ps_mean <- pseudomean(time = clindat$os_year, event = clindat$death_final, tmax = floor(max(clindat$os_year))) # negative values
