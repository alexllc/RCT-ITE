source("./bin/load_lib.R")


# trial_hte_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
#                 "NCT00460265", # CF takes a long time
#                 "NCT00041119_length", "NCT00041119_chemo",
#                 "NCT00003299", "NCT00119613")
trial_hte_ls <- c("NCT00113763")

min_mse_method <- read.csv("./res/crossfit_rloss/best_tau_estimators.csv")

for (trial in trial_hte_ls) {
    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]

    trial_best_method <- min_mse_method[min_mse_method$trial == trial, "best_tau_method"]
    tau <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", trial_best_method, "_tau_estimates.csv"))

    if (trial_best_method == "CF") {
        tau_vec <- tau$predictions
    } else {
        tau_vec <- tau[,1]
    }

    message("Performing cross-fitted XGBoost.")
    xgb_res <- cvboost(x = X, y = tau_vec, objective="reg:squarederror")
    saveRDS(xgb_res, file = paste0("./dat/xgb_model/", trial, "_xgb_model.rds"))
    message("XGBoost model saved.")
}