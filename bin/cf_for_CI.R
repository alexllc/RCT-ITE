source("./bin/load_lib.R")

min_mse_method <- read.csv("./res/best_tau_estimators.csv")
pos_report <- c("NCT00113763", "NCT00115765", "NCT00339183", "NCT00364013", "NCT00460265")
trial_choice <- filter(min_mse_method, trial %in% pos_report)
colnames(trial_choice)[1] <- "trialID"

for (j in 1:dim(trial_choice)[1]) {
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

    # Build causal forest
    cf <- causal_forest(X, Y, W, num.trees = 100000, tune.parameters = "all")
    cf_tau_pred <- predict(cf, estimate.variance =  TRUE)

    tau_df <- data.frame(cf_tau_pred, tau_est = tau[,1])

    # Extract standard error of predictions
    minid <- which(tau_df$tau_est == min(tau_df$tau_est))
    medid <- which.min(abs(tau_df$tau_est - median(tau_df$tau_est)))
    meanid <- which.min(abs(tau_df$tau_est - mean(tau_df$tau_est)))
    maxid <- which(tau_df$tau_est == max(tau_df$tau_est))

    tauci <- c(cf_tau_pred$variance.estimates[minid], cf_tau_pred$variance.estimates[medid], cf_tau_pred$variance.estimates[meanid], cf_tau_pred$variance.estimates[maxid])
    names(tauci) <- c("var.min", "var.med", "var.mean", "var.max")

    rtau <- c(summary(tau[,1])[c(1,3,4,6)], tauci, var(tau[,1]))
    
    
    if (j ==1){
        ci_df <- t(as.data.frame(rtau))
    } else {
        ci_df <- rbind(ci_df, t(as.data.frame(rtau)))
    }
}

ci_df <- cbind(select(trial_choice, c("trialID", "outcome_type")), ci_df)
report_ci <- data.table(ci_df) %>% mutate_if(is.numeric, ~round(.,2)) %>%  mutate(min_ci = paste0(Min., " (", Min. - 2*var.min, " to ", Min. + 2*var.min, ")")) %>% mutate(med_ci = paste0(Median, " (", Median - 2*var.med, " to ", Median + 2*var.med, ")")) %>% mutate(mean_ci = paste0(Mean, " (", Mean - 2*var.mean, " to ", Mean + 2*var.mean, ")")) %>% mutate(max_ci = paste0(Max., " (", Max. - 2*var.max, " to ", Max. + 2*var.max, ")")) %>% select(c("trialID", "outcome_type", "min_ci", "med_ci", "mean_ci", "max_ci"))

write.csv(report_ci, "./res/ite_tau_estimates/ITE_sum_CI.csv", row.names = FALSE)