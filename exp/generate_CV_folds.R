#' Script to generate the CV folds for each triral
#' Please reuse the same set of holdout vs training assignments for reproducible results

source("./bin/load_lib.R")
Q <- 4

# Extract trial names from list of scripts
trial_scripts <- list.files("./bin/load_RCT")
trial_ls <- trial_scripts[grep("^load*", trial_scripts)]
trial_ls <- unlist(lapply(strsplit(trial_ls, "_|\\."), function(X) X[2]))


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