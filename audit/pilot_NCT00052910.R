# Pilot study using the NCT00052910 dataset
# CALGB 80101: Phase III Intergroup Trial of Adjuvant Chemoradiation After Resection of Gastric or Gastroesophageal Adenocarcinoma
# RATIONALE: Drugs used in chemotherapy use different ways to stop tumor cells from dividing so they stop growing or die. Radiation therapy uses high-energy x-rays to damage tumor cells. Combining chemotherapy with radiation therapy after surgery may kill any remaining tumor cells following surgery. It is not yet known which chemotherapy and radiation therapy regimen is more effective in treating stomach or esophageal cancer.
# https://clinicaltrials.gov/ct2/show/study/NCT00052910

library(grf)
library(dplyr)
library(survival)
library(survminer)
library(timereg)
library(ggfortify) # plotting KM
library(pseudo)

df <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Multiple_Allianc_2002_213_NCT00052910/NCT00052910_D1_(EFFICACY).csv")
# univaraite KM anaylsis of arm I vs arm II
fit <- survfit(Surv(os_year, death_final) ~ TREAT_ASSIGNED, data = df)
autoplot(fit)

# multivariate cox PH 
mvmodel <- coxph(Surv(os_year, death_final) ~ T_stage + N_stage + M_stage + SEX_ID + agecat + Histologic_grade, data = df)  
 
# Check assumption of covariates *not* varying over time
aa_fit <- aareg(Surv(os_year, death_final) ~ T_stage + N_stage + M_stage + SEX_ID + agecat + Histologic_grade, data = df)

# use timereg to adjust for TVE
res <- comp.risk(Surv(os_year, death_final))

tx <- df$TREAT_ASSIGNED
outcome_col <- c("death_final", "os_year", "DFS_status", "DFS_year")
outcome <- dplyr::select(df, all_of(outcome_col))
covar <- dplyr::select(df, -c("MASK_ID", "TREAT_ASSIGNED", outcome_col))

cf <- causal_forest(
                    X = covar,
                    Y = comp_df$imputeyn_mean,
                    W = tx

)

pseudo_df <- data.frame(time = , event = )
cutoffs <- c(1, 3, 5, 10, floor(max(df$os_year)))
ps <- pseudosurv(time = df$os_year, event = df$death_final, tmax = cutoffs)
censored_ps <- cbind(df$death_final, df$os_year, ps$pseudo)

ps_mean <- pseudomean(time = df$os_year, event = df$death_final, tmax = floor(max(df$os_year)))
ext_psmean <- cbind(df$death_final, df$os_year, ps_mean)

#rearrange the data into a long data set
b <- NULL
for(it in 1:length(pseudo$time)){
b <- rbind(b,cbind(bmt,pseudo = pseudo$pseudo[,it],
tpseudo = pseudo$time[it],id=1:nrow(bmt)))
X = cbind(df$T_stage, df$N_stage)
X[is.na(X)] <- 0
yn <- imputeYn(X = X, Y = df$os_year, delta = df$death_final)


# Tail correction
efron_tail <- impute.survival(surv.time = df$os_year, censor = df$death_final)

# Comparison
OS <- df$death_final
time <- df$os_year
comp_df <- data.frame(OS, time, ps_mean, exp(efron_tail))

cor.test(comp_df$ps_mean, comp_df$exp.efron_tail.)