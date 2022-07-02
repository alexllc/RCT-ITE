# Script to exclusively find best R-loss parameters for building causal boosting models

source("./bin/load_lib.R")
source("./bin/hyperparameter_tuning.R")
source("./bin/mod_PTOforest.R")
source("./bin/mod_causalboost.R")

prob <- 0.5
Q <- 4


# Extract trial names from list of scripts
trial_scripts <- list.files("./bin/load_RCT")
trial_ls <- trial_scripts[grep("^load*", trial_scripts)]
trial_ls <- unlist(lapply(strsplit(trial_ls, "_|\\."), function(X) X[2]))

for (trial in trial_ls) {
    # For testing only:
    # Done 
    # NCT00339183, NCT00364013

    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
    
    print("Regular param testing.")
    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]

    outcome_list <- get(paste0(trial, "_outcomes"))

    for (outcome in outcome_list) {
        message(paste0(rep("=", 80)))
        message(paste0("Analyzing: ", outcome))
        message(paste0(rep("=", 80)))
        
        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }
        
        # We need to take out 1/Q sample as the holdout test data, the rest of the samples are used for training
        available_sample <- seq(1:dim(X)[1])
        holdout_sample <- sample(available_sample, floor(dim(X)[1] / Q))
        train_sample <- available_sample[!available_sample %in% holdout_sample]

        ho_X <- X[holdout_sample,]
        ho_W <- W[holdout_sample]
        ho_Y <- Y[holdout_sample]

        train_X <- X[train_sample,]
        train_W <- W[train_sample]
        train_Y <- Y[train_sample]

        cb_param <- find_cb_param(train_X, train_Y, train_W, num_search_rounds = 10)

        print(cb_param)

        write.csv(cb_param, paste0("./dat/cb_param/", trial, "_", imp_type_Y, "_", outcome, "_cb_param.csv"), row.names = TRUE)
    }
}
