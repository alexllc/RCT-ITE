source("./bin/load_lib.R")

linear_loss <- read.csv("./res_linear/compare_rloss.csv")
need_log <- filter(linear_loss, stack > 200)
trial_list <- unique(need_log$trial)

for (trial in trial_list) {

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

    outcome_list <- need_log$outcome[need_log$trial == trial]

    for (outcome in outcome_list) {

        message(paste0("Processing outcome: ", outcome))

        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efron+Yn"
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
        
        tau <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", imp_type_Y, "_", outcome, "_Rstack_tau_estimates.csv"))
        tau_vec <- tau[,1]

        rmv_id <- check_lm_rmv(X = X, Y = Y)

        X_rmvd <- X[-rmv_id,]
        Y_rmvd <- Y[-rmv_id]
        W_rmvd <- W[-rmv_id]

        message("Performing cross-fitted XGBoost.")
        xgb_res <- cvboost(x = X_rmvd, y = tau_vec, objective="reg:squarederror")
        saveRDS(xgb_res, file = paste0("./dat/xgb_model/", trial, "_", outcome, "_Rstack_xgb_model.rds"))
        message("XGBoost model saved.")
    }
}