# Script to exclusively find best R-loss parameters for building causal boosting models

source("./bin/load_lib.R")
source("./bin/hyperparameter_tuning.R")
source("./bin/mod_PTOforest.R")
source("./bin/mod_causalboost.R")

prob = 0.5
Q = 4
tuned_cb_param = TRUE
tuned_cm_param = TRUE
tune_ptof_param = TRUE
perform_xb = TRUE
# 
trial_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274")

for (trial in trial_ls) {
    # For testing only:
    # Done 
    # NCT00339183, NCT00364013

    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))

    source(paste0("./bin/load_CRC_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]


    # We need to take out 1/Q sample as the holdout test data, the rest of the samples are used for training
    holdout_sample <- sample(1:dim(X)[1], dim(X)[1] / Q)
    train_sample <- seq(1:dim(X)[1])[!(seq(1:dim(X)[1]) %in% holdout_sample)]
    ho_X <- X[holdout_sample,]
    ho_W <- W[holdout_sample]
    ho_Y <- Y[holdout_sample]

    train_X <- X[train_sample,]
    train_W <- W[train_sample]
    train_Y <- Y[train_sample]

    cb_param <- find_cb_param(train_X, train_Y, train_W, num_search_rounds = 10)

    print(cb_param)

    write.csv(cb_param, paste0("./dat/", trial, "_cb_param.csv"), row.names = TRUE)
    write.csv(holdout_sample, paste0("./dat/", trial, "_holdout_id.csv"), row.names = FALSE)
}
