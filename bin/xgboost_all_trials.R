source("./bin/load_lib.R")

# trial_hte_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
#                 "NCT00460265", # CF takes a long time
#                 "NCT00041119_length", "NCT00041119_chemo",
#                 "NCT00003299", "NCT00119613")
# trial_hte_ls <- c("NCT00113763")

min_mse_method <- read.csv("./res/best_tau_estimators.csv")
pos_report <- c("NCT00113763", "NCT00115765", "NCT00339183", "NCT00364013", "NCT00460265")
trial_choice <- filter(min_mse_method, trial %in% pos_report)
colnames(trial_choice)[1] <- "trialID"

# for (j in 1:dim(trial_choice)[1]) {

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
    
    trial_best_method <- trial_choice[j,2]

    tau <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", outcome, "_", trial_best_method, "_tau_estimates.csv"))

    tau_vec <- tau[,1]

    message("Performing cross-fitted XGBoost.")
    xgb_res <- cvboost(x = X, y = tau_vec, objective="reg:squarederror")
    saveRDS(xgb_res, file = paste0("./dat/xgb_model/", trial, "_", outcome, "_", trial_best_method, "_xgb_model.rds"))
    message("XGBoost model saved.")
# }