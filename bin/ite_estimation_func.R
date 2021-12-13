# Causal inference functions
library(causalToolbox) # X-learner
library(grf)
library(causalLearning) # Powers

library(magrittr)

avg_splits_res <- function(type_out = NULL, output_colnames = list()) {
    avg_df <- c()
    for (out_item in output_colnames) {
        j <- 1
        item_avg <- rowMeans(type_out[seq(j, splitRep * length(out_items), by = length(out_items))])
        j <- j + 1
        avg_df <-  cbind(avg_df, item_avg)
    } 
    avg_df <- as.data.frame(avg_df)
    colnames(avg_df) <- paste0("avg_", output_colnames)
    return(avg_df)
}

causal_X <- function(W = NULL, Y = NULL, X = NULL, splitRep = 10, ciBoostrap = 10, seed = NULL) {

    out_rf <- matrix(,nrow = length(W), ncol = 0)
    out_bart <- matrix(,nrow = length(W), ncol = 0)

    for (i in 1:splitRep) {
        message(paste0(">=====================>Repeat number: ", i))
        # assign group
        id <- 1:length(W)
        train_set <- sample(id, size = floor(length(W)/2))
        test_set <- id[!(id %in% train_set)]

        # train RF and BART using the first half
        xl_rf1 <- X_RF(feat = X[train_set,], tr = W[train_set], yobs = Y[train_set])
        xl_bart1 <- X_BART(feat = X[train_set,], tr = W[train_set], yobs = Y[train_set])

        # train RF and BART using the second half
        xl_rf2 <- X_RF(feat = X[test_set,], tr = W[test_set], yobs = Y[test_set])
        xl_bart2 <- X_BART(feat = X[test_set,], tr = W[test_set], yobs = Y[test_set])

        # Estimate and obtain confidence intervals for the two halves of RF and BART
        xl_rf1_ci <- CateCI(xl_rf1, X[test_set,], B = ciBoostrap)
        xl_rf2_ci <- CateCI(xl_rf2, X[train_set,], B = ciBoostrap)
        all_rf_ci <- data.frame(pred = double(), X.05 = double(), X.95 = double())
        all_rf_ci[test_set,] <- xl_rf1_ci
        all_rf_ci[train_set,] <- xl_rf2_ci

        xl_bart1_ci <- CateCI(xl_bart1, X[test_set,], B = ciBoostrap)
        xl_bart2_ci <- CateCI(xl_bart2, X[train_set,], B = ciBoostrap)
        all_bart_ci <- data.frame(pred = double(), X.05 = double(), X.95 = double())
        all_bart_ci[test_set,] <- xl_bart1_ci
        all_bart_ci[train_set,] <- xl_bart2_ci

        # save to split aggregating output matrix
        out_rf <- cbind(out_rf, all_rf_ci)
        out_bart <- cbind(out_bart, all_bart_ci)
    }

    # average the output over number of splits
    out_items <- c("pred", "X.05", "X.95")
    split_avg_rf <- avg_splits_res(out_rf, out_items)
    split_avg_bart <- avg_splits_res(out_bart, out_items)

    return(list(split_avg_rf, split_avg_bart))
}

causal_GRF <- function(W = NULL, Y = NULL, X = NULL, grf_W_hat = 0.5, grf_num_trees = 6000, split = FALSE, splitRep = 10, seed = NULL) {

    if (!split) { # estimating ITE without splitting data into train and test
        # splitting into test and train set is optional
        cf <- causal_forest(X = X, Y = Y, W = W, W.hat = grf_W_hat, num.trees = grf_num_trees)
        tau_hat <- predict(cf, estimate.variance = TRUE)
        ate_all <- average_treatment_effect(cf, target.sample = "all")
        tc <- test_calibration(cf)

        # add confidence interval for treatment effect predction
        tau_hat <- tau_hat %>% mutate(sigma = sqrt(variance.estimates),
                            ci_up = predictions + sigma*1.96,
                            ci_low = predictions - sigma*1.96,
                            sig_TE = ci_up > 0 & ci_low > 0 | ci_up < 0 & ci_low < 0)
    return(list(tau_hat, ate_all, tc))
    } else {

        id <- 1:length(W)
        out_tau_hat <- matrix(,nrow = length(W), ncol = 0)
        out_tc <- matrix(, nrow = 0, ncol = 8)
        out_ate <- matrix(, nrow = 0, ncol = 2)

        for(i in 1:splitRep) {
            train_set <- sample(id, size = floor(length(W)/2))
            test_set <- id[!(id %in% train_set)]


            cf1 <- causal_forest(X = X[train_set,], Y = Y[train_set], W = W[train_set], W.hat = grf_W_hat, num.trees = grf_num_trees)
            cf2 <- causal_forest(X = X[test_set,], Y = Y[test_set], W = W[test_set], W.hat = grf_W_hat, num.trees = grf_num_trees)

            # generate output
            tau_hat1 <- predict(cf1, newdata = X[test_set,], estimate.variance = TRUE)
            tau_hat2 <- predict(cf2, newdata = X[train_set,], estimate.variance = TRUE)
            all_tau_hat <- data.frame(predictions = double(), variance.estimates = double())
            all_tau_hat[train_set,] <- tau_hat1
            all_tau_hat[test_set,] <- tau_hat2

            all_tau_hat <- all_tau_hat %>% mutate(sigma = sqrt(variance.estimates),
                    ci_up = predictions + sigma*1.96,
                    ci_low = predictions - sigma*1.96,
                    sig_TE = ci_up > 0 & ci_low > 0 | ci_up < 0 & ci_low < 0)
            
            out_tau_hat <- cbind(out_tau_hat, all_tau_hat)

            # estimate ATE
            ate_all1 <- average_treatment_effect(cf1, target.sample = "all")
            ate_all2 <- average_treatment_effect(cf2, target.sample = "all")

            all_ate <- data.frame(estimate = double(), std_err = double())
            all_ate[1,] <- ate_all1
            all_ate[2,] <- ate_all2

            out_ate <- rbind(out_ate, all_ate)

            tc1 <- test_calibration(cf1)
            tc2 <- test_calibration(cf2)

            all_tc <- data.frame(mean_est = double(), diff_est = double(), mean_se = double(), diff_se = double(), mean_t = double(), diff_t =  double(), mean_p = double(), diff_p = double())

            all_tc[1,] <- c(tc1)
            all_tc[2,] <- c(tc2)

            out_tc <- rbind(out_tc, all_tc)

        }
        return(list(avg_splits_res(out_tau_hat, ), avg_splits_res(out_ate), avg_splits_res(out_tc)))
    }

}



#' fucntion to generate ITE estimation
#'
#' @param W binary treatment vector
#' @param Y outcome vector
#' @param X n * p covariate matrix

causality <- function(W = NULL, Y = NULL, X = NULL, seed = 123) {

    ## Estimate with X learner

}
