#' Script to generate the CV folds for each triral
#' Please reuse the same set of holdout vs training assignments for reproducible results

source("./bin/load_lib.R")
Q <- 4
# "NCT00003299", "NCT00041119_chemo","NCT00041119_length", "NCT00052910","NCT00079274", "NCT00113763","NCT00115765_oxa", "NCT00115765","NCT00119613", "NCT00339183_mt_KRAS", "NCT00339183_wt_KRAS", "NCT00364013", "NCT00460265" 
trial_ls <- c("NCT00364013_wtKRAS")

for (trial in trial_ls) {

    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
        X <- as.matrix(get(trial)[[1]])
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
        X <- as.matrix(get(trial)[[1]])
    }
    available_sample <- seq(1:dim(X)[1])
    ho_matrix <- matrix(nrow = floor(dim(X)[1] / Q), ncol = Q)
    for (iter in 1:Q) {
        holdout_sample <- sample(available_sample, floor(dim(X)[1] / Q))
        available_sample <- available_sample[!available_sample %in% holdout_sample]
        ho_matrix[,iter] <- holdout_sample
    }
    write.csv(ho_matrix, paste0("./dat/cv_ho_ids/cb_param_", trial, "_holdout_id.csv"), row.names = FALSE)

    # Generate for best tau estimators
    available_sample <- seq(1:dim(X)[1])
    ho_matrix <- matrix(nrow = floor(dim(X)[1] / Q), ncol = Q)
    for (iter in 1:Q) {
        holdout_sample <- sample(available_sample, floor(dim(X)[1] / Q))
        available_sample <- available_sample[!available_sample %in% holdout_sample]
        ho_matrix[,iter] <- holdout_sample
    }
    write.csv(ho_matrix, paste0("./dat/cv_ho_ids/est_tau_", trial, "_holdout_id.csv"), row.names = FALSE)

}