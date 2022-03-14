library(Counterfactual)
library(quantreg)

#Counterfactual distribution of X constructed by transformation of reference distribution
## Not run:
data(engel)
attach(engel)
counter_income <- mean(income)+0.75*(income-mean(income))
rqres <- counterfactual(foodexp~income, counterfactual_var=counter_income,nreg=100, transformation=TRUE, sepcore = TRUE, ncore=2)
## End(Not run)
# Wage decomposition: counterfactual and reference populations correspond to different groups
data(nlsw88)
attach(nlsw88)
lwage <- log(wage)
# method: logit
logitres<-counterfactual(lwage~tenure+ttl_exp+grade, group=union, treatment=TRUE, decomposition=TRUE, method="logit", noboot=FALSE, sepcore = TRUE,ncore=2)



source("/home/alex/Documents/lab/RCT-ITE/audit/model_eval/load_lib.R")
source("/home/alex/Documents/lab/RCT-ITE/audit/model_eval/load_data.R")


efron_tail <- exp(impute.survival(surv.time = test_clin_df$PFS_mo, censor = test_clin_df$PFSCR))
pto_model <- PTOforest(x = X, tx = W, y = efron_tail, num.trees = 10000)
pto_tau <- predict(pto_model, X)

attach(X)

rctres <- counterfactual(efron_tail~LIVERMET+AGE+SEX+RACE+B_ECOG+DIAGTYPE, group = W, method = "logit", nreg=100, transformation=TRUE,  decomposition = TRUE, sepcore = TRUE, ncore=4)