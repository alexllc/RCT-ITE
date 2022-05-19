source("./bin/load_lib.R")

min_mse_method <- read.csv("./res/best_tau_estimators.csv")
trial_choice <- filter(min_mse_method, mse < 150)
colnames(trial_choice)[1] <- "trialID"

cstat_res <- data.frame(trial = character(), outcome = character(), cstat = numeric())
for (j in 1:dim(min_mse_method)[1]) {
    trial <- min_mse_method[j,1]

    message(paste0(rep("=", 80)))
    message(paste0("Loading: ", trial))
    message(paste0(rep("=", 80)))
    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    message("Data loaded.")
    
    outcome <- min_mse_method[j,4]
    tau_method <- min_mse_method[j,2]

    # load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]
    Y_list <- get(paste0(outcome, "_Y_list"))

    imp_type_Y <- "efronYn"
    Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

    if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
        Y <- as.numeric(Y_list[[1]])
        imp_type_Y <- "efron"
    }

    tau <- try(read.csv(paste0("./res/ite_tau_estimates/", trial, "_", outcome, "_", tau_method , "_tau_estimates.csv")))
    if (class(tau) == "try-error")
        next
    te_df <- data.frame(id = sel_SUBJID, tau = tau, tx = W, obs_Y = Y)
    te_df$tx <- ifelse(te_df$tx == 0, "control", "intervention")
    colnames(te_df) <- c("subj", "te", "arm", "obs_Y")
    c_df <- spread(te_df, arm, te)

    cctrl <- filter(c_df, !is.na(control))
    cctrl <- cctrl[order(cctrl$control),]
    cctrl$intervention <- NULL

    ctrt <- filter(c_df, !is.na(intervention))
    ctrt <- ctrt[order(ctrt$intervention),]
    ctrt$control <- NULL

    compare_length <- ifelse(dim(cctrl)[1] > dim(ctrt)[1], dim(ctrt)[1], dim(cctrl)[1])

    if (dim(cctrl)[1] != compare_length) {
        random_pick <- sample(1:dim(cctrl)[1], compare_length)
        cctrl <- cctrl[random_pick,]
        cctrl <- cctrl[order(cctrl$control),]
    } else {
        random_pick <- sample(1:dim(ctrt)[1], compare_length)
        ctrt <- ctrt[random_pick,]
        ctrt <- ctrt[order(ctrt$intervention),]
    }

    cstat_df <- cbind(ctrt, cctrl)
    colnames(cstat_df) <- c("trt_id", "trt_pred_te", "trt_obs", "ctrl_id", "ctrl_pred_te", "ctrl_obs")
    cstat_df <- data.table(cstat_df)
    cstat_df <- cstat_df %>% rowwise() %>% mutate(pred_diff = mean(trt_pred_te, ctrl_pred_te)) %>% mutate(obs_diff = trt_obs - ctrl_obs)

    concordance <- 0
    total <- 0
    for (i in 1:compare_length) {
        message(paste0("Iteration: ", i, " of ", compare_length))
        for (j in (1:compare_length)[1:compare_length != i]) {
            if ( ((cstat_df[i,7] < cstat_df[j,7]) & (cstat_df[i,8] < cstat_df[j,8])) | ((cstat_df[i,7] > cstat_df[j,7]) & (cstat_df[i,8] > cstat_df[j,8])))
                concordance <- concordance + 1
            total <- total + 1
        }
    }
    cstat <- concordance / total
    cstat_res <- rbind(cstat_res, c(trial, outcome, cstat))
}
    # NCT00052910 CBoost OS
    # 0.5045822

    # NCT00339183 stack PFS
    # 0.5314068
write.csv(cstat_res, "cstat.csv", row.names = FALSE)