#' Script to generate the CV folds for each triral
#' Please reuse the same set of holdout vs training assignments for reproducible results

Q <- 4

trial_ls <- c("NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274", "NCT00460265", "NCT00041119",  "NCT00003299", "NCT00119613", "NCT00041119", "NCT00364013_KRASe2", "NCT00364013_biom")

for (trial in trial_ls) {

    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
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