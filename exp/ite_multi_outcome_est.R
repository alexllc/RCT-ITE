# Script to find hetereogeneity based on the best selected tau values
source("./bin/load_lib.R")
source("./bin/rstack.R")

Q <- 4

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



# Extract trial names from list of scripts
trial_scripts <- list.files("./bin/load_RCT")
trial_ls <- trial_scripts[grep("^load*", trial_scripts)]
trial_ls <- unlist(lapply(strsplit(trial_ls, "_|\\."), function(X) X[2]))

# load linear R-loss values
linear_loss <- read.csv("./res_linear/compare_rloss.csv")
need_log <- filter(linear_loss, stack < 200) # even low R-loss trials will now be treated with log transformation
trial_log <- unique(need_log$trial)
# trial_log <- trial_log[5:length(trial_log)]

    for (trial in trial_log) {
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
        outcome_list <- get(paste0(trial, "_outcomes"))

        # only run outcomes that have higher than 200 R-losses
        outcome_list <- need_log$outcome[need_log$trial == trial]

        for (outcome in outcome_list) {

            message(paste0("Processing outcome: ", outcome))

            Y_list <- get(paste0(outcome, "_Y_list"))

            imp_type_Y <- "efronYn"
            Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

            if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
                Y <- as.numeric(Y_list[[1]])
                imp_type_Y <- "efron"
            }

            if(Y != "RSP") {
                Y <- log(Y)
                if (any(is.infinite(Y)))
                    Y[which(is.infinite(Y))] <- 0
            }
            
            rmv_id <- check_lm_rmv(X = X, Y = Y)

            X_rmvd <- X[-rmv_id,]
            Y_rmvd <- Y[-rmv_id]
            W_rmvd <- W[-rmv_id]

            tau <- perform_rstack(Y = Y_rmvd, X = X_rmvd, W = W_rmvd, trial_name = trial, outcome = outcome, imp_type_Y = imp_type_Y, tuned_cb_param = FALSE, tuned_cm_param = TRUE, tune_ptof_param = TRUE, perform_xb = FALSE, perform_cb = FALSE)
            write.csv(tau, paste0("./res/ite_tau_estimates/", trial, "_", imp_type_Y, "_", outcome, "_Rstack_tau_estimates.csv"), row.names = FALSE)
        } # outcome type
        
    } # trials

