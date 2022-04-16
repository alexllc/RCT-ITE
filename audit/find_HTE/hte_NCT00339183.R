source("./bin/load_lib.R")
source("./bin/hyperparameter_tuning.R")
source("./bin/mod_PTOforest.R")
source("./bin/mod_causalboost.R")


trial <- "NCT00339183"

source(paste0("./bin/load_CRC_RCT/load_", trial, ".R"))
X <- as.matrix(get(trial)[[1]])
Y <- get(trial)[[2]][[1]]
W <- get(trial)[[3]]

# CF has lowest MSE

cf <- causal_forest(X, Y, Y, num.trees = 100000, tune.parameters = "all")
cf_tau_pred <- predict(cf, estimate.variance =  TRUE)
while (var(cf_tau_pred$predictions) == 0) { # tuned forests sometimes give equal predictions, if that's the case, rebuild forest
    cf <- causal_forest(X, Y, W, num.trees = 100000, tune.num.trees = 100000, tune.parameters = "all")
    cf_tau_pred <- predict(cf, estimate.variance =  TRUE)
}


# Find if HTE is significant

test_hte <- data.frame(Y  = Y, Z = W)

hte_no_adj <- detect_idiosyncratic(formula = Y ~ Z, data = test_hte, plugin = TRUE, test.stat = "SKS.stat",  tau.hat = cf_tau_pred, n.cores = 40)

# FRT Plug-in Test for 
# Treatment Effect Heterogeneity 
#   Statistic P-Value (Sweep) P-Value (Plug-In)
# 1 0.1833274          1.0001            0.1381