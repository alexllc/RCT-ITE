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
    trial_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
                    "NCT00460265",
                    "NCT00041119_length", "NCT00041119_chemo",
                    "NCT00003299", "NCT00119613")

    for (trial in trial_ls) {

        for (iter in 1:best_tau_fold) {
            r_loss <- read.csv(paste0("./res/crossfit_rloss/", trial, "_model_sel_", iter, "_res.csv"))
            if (iter == 1){
                sum_tbl <- data.frame(method_name = r_loss$X)
                sum_tbl <- cbind(sum_tbl, r_loss$mse)
            } else {
                sum_tbl <- cbind(sum_tbl, r_loss$mse)
            }
        }
        trial_cv_mse <- rowMeans(sum_tbl[,2:(best_tau_fold+1)])

        if (trial == "NCT00364013") { # create new data frame for the first trial
            mean_rloss <- data.frame(method_name = sum_tbl$method_name)
            mean_rloss <- cbind(mean_rloss, trial_cv_mse)
            colnames(mean_rloss)[which(colnames(mean_rloss) == "trial_cv_mse")] <- trial
        } else {
            mean_rloss <- cbind(mean_rloss, trial_cv_mse)
            colnames(mean_rloss)[which(colnames(mean_rloss) == "trial_cv_mse")] <- trial
        }
    }

    min_mse_method <- data.frame()
    for (col in 2:dim(mean_rloss)[2]) {
        min_mse_method <- rbind(min_mse_method, c(colnames(mean_rloss)[col], mean_rloss$method_name[which.min(mean_rloss[,col])], min(mean_rloss[,col], na.rm = TRUE)))
    }

    colnames(min_mse_method) <- c("trial", "best_tau_method", "mse")

    write.csv(min_mse_method, "best_tau_estimators.csv", row.names = FALSE)
} else {
    min_mse_method <- read.csv("./res/crossfit_rloss/best_tau_estimators.csv")
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
trial_hte_ls <- c("NCT00113763")

hte_df <- data.frame(trial = character(), stat = numeric(), pval_sweep = numeric(), pval_plugin = numeric())

for (trial in trial_hte_ls) {

    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
    
    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]

    if (trial == "NCT00364013"){ # rboosting

        tau <- perform_rlearn(X = X, W = W, Y = Y, tauMethod = "rboost")

    } else if (trial == "NCT00079274") { # causal forest

        cf <- causal_forest(X, Y, W, num.trees = 100000, tune.parameters = "all")
        tau <- predict(cf, estimate.variance =  TRUE)
        while (var(tau$predictions) == 0) { # tuned forests sometimes give equal predictions, if that's the case, rebuild forest
            cf <- causal_forest(X, Y, W, num.trees = 100000, tune.parameters = "all")
            tau <- predict(cf, estimate.variance =  TRUE)
        }

    } else if (trial == "NCT00041119_length") { # causal MARS

        cm_cv_performance <- find_cm_param(x = X, tx = W, y = Y, verbose = TRUE)
        cm_param_used <- as.list(cm_cv_performance[1,1:3]) # choose smallest first
        cm <- try(do.call(causalMARS, append(list(x = X, tx = W, y = Y), cm_param_used)))
        
        tau <- try(predict(cm, X))


    } else if (trial %in% c("NCT00113763", "NCT00460265")) { # r-stack
        tau <- perform_rstack(Y = Y, X = X, W = W, trial_name = trial)

    } else { # causal boosting for the rest
        cb_param_tune_res <- read.csv(paste0("./dat/cb_param/", trial, "_cb_param.csv")) # if skipping parameters tuning, we will use the previously tuned parameter for this trial
        cb_param_tune_res <- filter(cb_param_tune_res, num_trees_min_effect != 1) # avoid using parameteres with only 1 boosting tree
        cb_param <- cb_param_tune_res[which.min(cb_param_tune_res$mean_cvm_effect),]
        cb_param <- list(num.trees = cb_param$num_trees_min_effect, maxleaves = cb_param$max_leaves, eps = cb_param$eps, splitSpread = cb_param$split_spread)

        cb <- do.call(causalBoosting, append(list(x = X, tx = W, y = Y), as.list(cb_param)))
        cb_tau_pred <- predict(cb, newx = X)
        tau <- cb_tau_pred[,cb_param$num.trees] # choose the prediction given by the maximum number of trees used 

    }

    write.csv(tau, paste0(trial, "_", min_mse_method[min_mse_method$trial == trial, "best_tau_method"], "_tau_estimates.csv"), row.names = FALSE)
    
    # hte <- detect_idiosyncratic(formula = Y ~ W, data = data.frame(Y, W), plugin = TRUE, tau.hat = tau, test.stat = "SKS.stat")

    # print(hte)
    
    # hte_df <- rbind(hte_df, c(trial, hte$statistic, hte$p.value, hte$p.value.plug))
    # write.csv(hte_df, paste0(trial, "_tmp_hte.csv"), row.names = FALSE)
    
}
# colnames(hte_df)  <- c("trial", "stat", "pval_sweep", "pval_plugin")

# write.csv(hte_df, "./res/HTE_results.csv", row.names = FALSE)