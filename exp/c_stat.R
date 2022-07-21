source("./bin/load_lib.R")

calcualte_c_benefit <- function(Y = NULL, tau = NULL, W = NULL) {
    te_df <- data.frame(id = 1:length(Y), tau = tau, tx = W, obs_Y = Y)
    te_df$tx <- ifelse(te_df$tx == 0, "control", "intervention")
    colnames(te_df) <- c("subj", "te", "arm", "obs_Y")
    c_df <- spread(te_df, arm, te)

    cctrl <- filter(c_df, !is.na(control))
    cctrl <- cctrl[order(cctrl$control),]
    cctrl$intervention <- NULL
    cctrl$control <- round(cctrl$control, digits = 0)

    ctrt <- filter(c_df, !is.na(intervention))
    ctrt <- ctrt[order(ctrt$intervention),]
    ctrt$control <- NULL
    ctrt$intervention <- round(ctrt$intervention, digits = 0)

    ctrt <- cbind(ctrt, cctrl[match(ctrt$intervention, cctrl$control),])
    ctrt <- ctrt[complete.cases(ctrt),]

    colnames(ctrt) <- c("trt_id", "trt_obs", "trt_pred_te", "ctrl_id", "ctrl_obs", "ctrl_pred_te")
    cstat_df <- data.table(ctrt)
    cstat_df <- cstat_df %>% rowwise() %>% mutate(mean_te = mean(trt_pred_te, ctrl_pred_te)) %>% mutate(obs_diff = trt_obs - ctrl_obs)

    c_benefit <- concordancefit(cstat_df$obs_diff, cstat_df$mean_te)
    return(c_benefit)
}



Q <- 4

# Extract trial names from list of scripts
ite_res_files <- list.files("./res/ite_tau_estimates/backup/")
trial_ls <- unique(unlist(lapply(strsplit(ite_res_files, "_|\\."), function(X) X[1])))

cstat_list <- list()

for (trial in trial_ls) {
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
    
    trial_files <- ite_res_files[grep(paste0("^", trial,"*"), ite_res_files)]
    outcome_list <- unlist(lapply(strsplit(trial_files, "_|\\."), function(X) X[3]))

    for (outcome in outcome_list) {

        message(paste0("Processing outcome: ", outcome))

        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }
        Y_lasso <- cv.glmnet(X, Y, keep = TRUE, family = "gaussian")
        Y_hat_lasso <- predict(Y_lasso, newx = X)[,1]

        tau_pred <- read.csv(paste0("./res/ite_tau_estimates/backup/", trial, "_", imp_type_Y, "_", outcome, "_Rstack_tau_estimates.csv"))

        tau_vec <- tau_pred[,1]

        cstat <- calcualte_c_benefit(Y = (Y - Y_hat_lasso)[1:length(tau_vec)], tau = tau_vec, W = W[1:length(tau_vec)])
        cstat$count <- NULL
        cstat_list <- append(cstat_list, list(c(trial, outcome, as.list(cstat))))
        message("C for benefit is: ")
        print(cstat)
       
        }
}
    # NCT00052910 CBoost OS
    # 0.5045822

    # NCT00339183 stack PFS
    # 0.5314068

cstat_sum <- do.call(rbind.data.frame, cstat_list)

colnames(cstat_sum) <- c("trial", "outcome", "Cbenefit", "n", "var", "cvar")
write.csv(cstat_sum, "./res/counter_factual_cstat.csv", row.names = FALSE)


sig_trial_outcome$key <- paste0(sig_trial_outcome$trial, sig_trial_outcome$V2)

cstat_sum$key <- paste0(cstat_sum$trial, cstat_sum$outcome)

sig_trial_outcome <- left_join(sig_trial_outcome, cstat_sum, by = c("key"))

non_sig_cstat <- cstat_sum[!(cstat_sum$key %in% sig_trial_outcome$key),]

write.csv(sig_trial_outcome, "./res/sig_trial_cstat.csv", row.names = FALSE)
write.csv(non_sig_cstat, "./res/non_sig_trial_cstat.csv", row.names = FALSE)