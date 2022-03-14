library(gamlss)
library(dplyr)
library(survival)
library(survminer)

library(causalLearning) # Powers

data(abdom)
mod<-gamlss(y~pb(x),sigma.fo=~pb(x),family=BCT, data=abdom, method=mixed(1,20))
plot(mod)
rm(mod)


clin_df <- read.csv("./dat/NCT00364013_adls.csv")
biom <- read.csv("./dat/NCT00364013_biom.csv")

tested_subj <- Reduce(union, lapply(biom[,c(2:dim(biom)[2])], function(X) which(is.na(X))))
test_clin_df <- clin_df[clin_df$SUBJID %in% tested_subj,]

dim(test_clin_df)

test_clin_df$PFS_mo <- test_clin_df$PFSDYCR / 30.4167
test_clin_df$OS_mo <- test_clin_df$DTHDY / 30.4167

X <- select(test_clin_df, all_of(c("LIVERMET", "DIAGMONS", "AGE",  "SEX", "B_WEIGHT", "B_HEIGHT", "RACE",  "B_ECOG", "HISSUBTY", "B_METANM", "DIAGTYPE")))
W <- test_clin_df$ATRT
source("./audit/original_surv.r")
efron_tail <- exp(impute.survival(surv.time = test_clin_df$PFS_mo, censor = test_clin_df$PFSCR)) # no need to convert to log scale

pto_model <- PTOforest(x = X, tx = W, y = efron_tail, num.trees = 10000)
pto_tau <- predict(pto_model, X)

mod <- gamlss(pto_tau~LIVERMET+AGE+SEX+B_WEIGHT+B_WEIGHT+RACE+B_ECOG+HISSUBTY+B_METANM, data = X, sigma.formula = ~LIVERMET+AGE+SEX+B_WEIGHT+B_WEIGHT+RACE+B_ECOG+HISSUBTY+B_METANM, family = BCT, method=mixed(1,20))
plot(mod)