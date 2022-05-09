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
# CRC:"NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274"
# HNC: "NCT00460265"
# BC: "NCT00041119"
# SCLC: "NCT00003299", "NCT00119613", "NCT00041119"

trial_ls <- c("NCT00364013_KRASe2", "NCT00364013_biom")

for (trial in trial_ls) {
    # For testing only:
    # Done 
    # NCT00339183, NCT00364013

    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))

    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    print("Regular param testing.")
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]

    outcome_list <- get(paste0(trial, "_outcomes"))

    for (outcome in outcome_list) {
        Y_list <- get(paste0(outcome, "_Y_list"))

        imputation_methods <- c("efron", "efron+Yn", "pseudo")
        for (imp_type_Y in 1:3) {
    
            if (is.na(imp_type_Y)) {
                message(paste0("Imputation type ", imputation_methods[imp_type_Y], " is not available, skipping."))
                next
            } else {
                Y <- Y_list[[imp_type_Y]]
            }

            ho_matrix <- read.csv(paste0("./dat/cv_ho_ids/cb_param_", trial, "_holdout_id.csv"))
            for (iter in 1:Q) {

                # We need to take out 1/Q sample as the holdout test data, the rest of the samples are used for training
                holdout_sample <- ho_matrix[,iter]
                available_sample <- seq(1:dim(X)[1])
                train_sample <- available_sample[!available_sample %in% holdout_sample]

                ho_X <- X[holdout_sample,]
                ho_W <- W[holdout_sample]
                ho_Y <- Y[holdout_sample]

                train_X <- X[train_sample,]
                train_W <- W[train_sample]
                train_Y <- Y[train_sample]

                cb_param <- find_cb_param(train_X, train_Y, train_W, num_search_rounds = 10)

                print(cb_param)

                write.csv(cb_param, paste0("./dat/cb_param/", trial, imputation_methods[imp_type_Y], "_", outcome, "_", iter, "_cb_param.csv"), row.names = TRUE)
            }
        }
    }
}
