# Script to find hetereogeneity based on the best selected tau values
source("./bin/load_lib.R")
source("./bin/rstack.R")

import_mse_res <- FALSE

# Function to perform r-learning over Q fold cross-fitting
perform_rlearn <- function(X, W, Y, Q = 4, tauMethod = NULL, prob = 0.5) {
    
    available_sample <- 1:dim(X)[1]
    cf_tau <- data.frame(ID = available_sample)
    rownames(cf_tau) <- cf_tau$ID

    for (iter in 1:Q) {
        
        if (iter != Q) {
            ho_id <- sample(available_sample, floor(dim(X)[1] / Q))
            available_sample <- available_sample[!(available_sample %in% ho_id)]
        } else {
            ho_id <- available_sample
        }

        train_id <- seq(1:dim(X)[1])[!(seq(1:dim(X)[1]) %in% ho_id)]

        if (tauMethod == "rboost") {
            mhat <- get_mhat(X[train_id,], W[train_id], Y[train_id], X[ho_id,], Y[ho_id])
            model <- do.call(tauMethod, list(x = X[train_id,], w = W[train_id], y = Y[train_id], p_hat = prob, m_hat = mhat, nthread = 40))
            ho_tau_pred <- predict(model, X[ho_id,])
        } else if (tauMethod == "causal_forest") {
            model <- causal_forest(X = X[train_id,], Y = Y[train_id], W[train_id], num.trees = 100000, tune.parameters = "all")
            ho_tau_pred <- predict(model, newdata = X[ho_id,], estimate.variance =  TRUE)
            
            # tuned forests sometimes give equal predictions, if that's the case, rebuild forest
            while (var(ho_tau_pred$predictions) == 0) { 
                model <- causal_forest(X = X[train_id,], Y = Y[train_id], W[train_id], num.trees = 100000, tune.parameters = "all")
                ho_tau_pred <- predict(model, newdata = X[ho_id,], estimate.variance =  TRUE)
            }
        }
        cf_tau[ho_id,] <- ho_tau_pred # save out of bag tau estimation to dataframe
    }
    return(cf_tau[,1])
}

perform_ptof <- function(X, W, Y, outcome = NULL, Q = 4, prob = 0.5) {

    available_sample <- 1:dim(X)[1]
    ptof_tau <- data.frame(id = available_sample)
    rownames(ptof_tau) <- ptof_tau$ID

    for (iter in 1:Q) {

        if (iter != Q) {
            ho_id <- sample(available_sample, floor(dim(X)[1] / Q))
            available_sample <- available_sample[!(available_sample %in% ho_id)]
        } else {
            ho_id <- available_sample
        }
        train_id <- seq(1:dim(X)[1])[!(seq(1:dim(X)[1]) %in% ho_id)]

        ho_X <- X[ho_id,]
        ho_W <- W[ho_id]
        ho_Y <- Y[ho_id]

        train_X <- X[train_id,]
        train_W <- W[train_id]
        train_Y <- Y[train_id]

        ptof_param_search <- find_ptof_param(x = train_X, y = train_Y, tx = train_W, validation_fold = 4, num_search_rounds = 50)
        ptof_param <- ptof_param_search[1, 1:3]
        ptof_param <- list(num.trees = ptof_param$num_trees, mtry = ptof_param$mtry, min.node.size = ptof_param$min_node_size)
        ptof <- try(do.call(PTOforest, append(list(x = train_X, y = train_Y, tx = train_W, postprocess = FALSE, verbose = TRUE), ptof_param))) # sometimes PTO forest can fail to produce some estimates with holdout data
        write.csv(ptof_param_search, paste0("./res/params/", trial, "_", outcome, "_", iter,  "_ptof_params_USED.csv"), row.names = FALSE)
        ptof_tau_pred <- try(predict(ptof, ho_X))
        ptof_tau[ho_id,] <- ptof_tau_pred # save out of bag tau estimation to dataframe
    }
    return (ptof_tau[,1])
}

perform_cb <- function(X, W, Y, outcome = NULL, Q = 4, prob = 0.5) {
    available_sample <- 1:dim(X)[1]
    ptof_tau <- data.frame(id = available_sample)
    rownames(ptof_tau) <- ptof_tau$ID

    for (iter in 1:Q) {

        if (iter != Q) {
            ho_id <- sample(available_sample, floor(dim(X)[1] / Q))
            available_sample <- available_sample[!(available_sample %in% ho_id)]
        } else {
            ho_id <- available_sample
        }
        train_id <- seq(1:dim(X)[1])[!(seq(1:dim(X)[1]) %in% ho_id)]

        ho_X <- X[ho_id,]
        ho_W <- W[ho_id]
        ho_Y <- Y[ho_id]

        train_X <- X[train_id,]
        train_W <- W[train_id]
        train_Y <- Y[train_id]

        ptof_param_search <- find_ptof_param(x = train_X, y = train_Y, tx = train_W, validation_fold = 4, num_search_rounds = 50)
        ptof_param <- ptof_param_search[1, 1:3]
        ptof_param <- list(num.trees = ptof_param$num_trees, mtry = ptof_param$mtry, min.node.size = ptof_param$min_node_size)
        ptof <- try(do.call(PTOforest, append(list(x = train_X, y = train_Y, tx = train_W, postprocess = FALSE, verbose = TRUE), ptof_param))) # sometimes PTO forest can fail to produce some estimates with holdout data
        write.csv(ptof_param_search, paste0("./res/params/", trial, "_", outcome, "_", iter,  "_ptof_params_USED.csv"), row.names = FALSE)
        ptof_tau_pred <- try(predict(ptof, ho_X))
        ptof_tau[ho_id,] <- ptof_tau_pred # save out of bag tau estimation to dataframe
    }
    return (ptof_tau[,1])
}

# Function to determine whether to use boosting or LASSO to derive mhat
get_mhat <- function(X_train, W_train, Y_train, X_ho, Y_ho) {

    Y_boost <- cvboost(X_train, Y_train, objective = "reg:squarederror", nthread = 40) # non-binary outcome
    Y_hat_boost <- predict(Y_boost, newx = X_ho)

    Y_lasso <- cv.glmnet(X_train, Y_train, keep = TRUE, family = "gaussian")
    Y_hat_lasso <- predict(Y_lasso, newx = X_ho)[,1]

    # which method has a smaller CV error?
    print("RMSE of m hat estimated by boosting vs LASSO:")
    print(round(c(RMSE(Y_hat_boost, Y_ho), RMSE(Y_hat_lasso, Y_ho)), 4))
    
    if (RMSE(Y_hat_boost, Y_ho) < RMSE(Y_hat_lasso, Y_ho)) {
        return(Y_hat_boost)
    } else {
        return(Y_hat_lasso)
    }
}

if (import_mse_res) {

best_tau_fold <- 4
trial_ls <- c("NCT00052910", "NCT00113763", "NCT00460265","NCT00364013",  "NCT00115765", "NCT00339183")

    for (trial in trial_ls) {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
        outcome_list <- get(paste0(trial, "_outcomes"))

        for (outcome in outcome_list) {
            imp_type_Y <- "efron+Yn"
            Y_list <- get(paste0(outcome, "_Y_list"))
            Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

            if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
                Y <- as.numeric(Y_list[[1]])
                imp_type_Y <- "efron"
            }

            for (iter in 1:best_tau_fold) {
                r_loss <- read.csv(paste0("./res/crossfit_rloss/", trial, "_model_sel_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"))
                if (iter == 1){
                    sum_tbl <- data.frame(method_name = r_loss$X)
                    sum_tbl <- cbind(sum_tbl, r_loss$mse)
                } else {
                    sum_tbl <- cbind(sum_tbl, r_loss$mse)
                }
            }
            trial_cv_mse <- rowMeans(sum_tbl[,2:(best_tau_fold+1)])

            if (trial == "NCT00052910" & outcome == "OS") { # create new data frame for the first trial
                mean_rloss <- data.frame(method_name = sum_tbl$method_name)
                mean_rloss <- cbind(mean_rloss, trial_cv_mse)
                colnames(mean_rloss)[which(colnames(mean_rloss) == "trial_cv_mse")] <- paste0(trial, "_", outcome)
            } else {
                mean_rloss <- cbind(mean_rloss, trial_cv_mse)
                colnames(mean_rloss)[which(colnames(mean_rloss) == "trial_cv_mse")] <- paste0(trial, "_", outcome)
            }
        } # outcome type
    } # trials

    min_mse_method <- data.frame()
    for (col in 2:dim(mean_rloss)[2]) {
        min_mse_method <- rbind(min_mse_method, c(colnames(mean_rloss)[col], mean_rloss$method_name[which.min(mean_rloss[,col])], min(mean_rloss[,col], na.rm = TRUE)))
    }
    min_mse_method$outcome_type <- unlist(lapply(min_mse_method$trial, function(X) strsplit(X, "_")[[1]][2]))
    min_mse_method$trial <- unlist(lapply(min_mse_method$trial, function(X) strsplit(X, "_")[[1]][1]))
    min_mse_method$mse <- as.numeric(min_mse_method$mse)
    colnames(min_mse_method) <- c("trial", "best_tau_method", "mse", "outcome_type")

    write.csv(min_mse_method, "./res/best_tau_estimators.csv", row.names = FALSE)
} else {
    min_mse_method <- read.csv("./res/best_tau_estimators.csv")
}
#
# Calculate HTE using the best tau estimator 
#
# Provide propensity score for RCT
prob <- 0.5

# trial_hte_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
#                 "NCT00460265", # CF takes a long time
#                 "NCT00041119_length", "NCT00041119_chemo",
#                 "NCT00003299", "NCT00119613")

trial_choice <- filter(min_mse_method, mse < 150)
colnames(trial_choice)[1] <- "trialID"



for (j in 5:dim(trial_choice)[1]) {

    trial <- trial_choice[j,1]
    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
    
    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]

    outcome <- trial_choice[j,4]

    message(paste0("Processing outcome: ", outcome))

    Y_list <- get(paste0(outcome, "_Y_list"))

    imp_type_Y <- "efron+Yn"
    Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

    if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
        Y <- as.numeric(Y_list[[1]])
        imp_type_Y <- "efron"
    }
    
    best_ite_model <- trial_choice[j,2]

    if (best_ite_model == "PTOForest"){ # PTOForest
        message("Building PTOForest.")
        tau <- perform_ptof(X, W, Y, outcome = outcome)

    } else if (best_ite_model == "CBoost" ) { # causal forest
        message("Building CBoosting")
        tau <- perform_cb(Y = Y, X = X, W = W)

    } else { # stack for the rest
        message("Performing RStack")
        tau <- perform_rstack(Y = Y, X = X, W = W, trial_name = trial, outcome = outcome, imp_type_Y = imp_type_Y)
    }

    write.csv(tau, paste0(trial, "_", outcome, "_", best_ite_model, "_tau_estimates.csv"), row.names = FALSE)
}

