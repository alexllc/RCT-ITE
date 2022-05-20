source("./bin/load_lib.R")

calcualte_c_benefit <- function(Y = NULL, tau = NULL, W = NULL) {
    te_df <- data.frame(id = 1:length(Y), tau = tau, tx = W, obs_Y = Y)
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

    cstat_df <- as.data.frame(cstat_df[,7:8])

    concordance <- 0
    total <- 0
    for (p in 1:compare_length) {
        # message(paste0("Iteration: ", p, " of ", compare_length))
        for (q in (1:compare_length)[1:compare_length != p]) {
            if ( ((cstat_df[p,1] < cstat_df[q,1]) & (cstat_df[p,2] < cstat_df[q,2])) | ((cstat_df[p,1] > cstat_df[q,1]) & (cstat_df[p,2] > cstat_df[q,2])))
                concordance <- concordance + 1
            total <- total + 1
        }
    }
    return(cstat <- concordance / total)
}


min_mse_method <- read.csv("./res/best_tau_estimators.csv")
pos_report <- c("NCT00113763", "NCT00115765", "NCT00339183", "NCT00364013", "NCT00460265")
trial_choice <- filter(min_mse_method, mse < 150 & trial %in% pos_report)
colnames(trial_choice)[1] <- "trialID"

cstat_res <- data.frame(trial = character(), outcome = character(), cstat = numeric())
for (j in 1:dim(trial_choice)[1]) {
    trial <- trial_choice[j,1]

    message(paste0(rep("=", 80)))
    message(paste0("Loading: ", trial))
    message(paste0(rep("=", 80)))
    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    message("Data loaded.")
    
    outcome <- trial_choice[j,4]
    tau_method <- trial_choice[j,2]

    message(paste0("Outcome: ", outcome))

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

    tau_vec <- tau[,1]

    cstat <- calcualte_c_benefit(Y = Y, tau = tau_vec, W = W)
    cstat_res <- rbind(cstat_res, c(trial, outcome, cstat))

    message("C for benefit is: ")
    print(cstat)
    print(cstat_res)
    
}
    # NCT00052910 CBoost OS
    # 0.5045822

    # NCT00339183 stack PFS
    # 0.5314068

colnames(cstat_res) <- c("trial", "outcome", "Cbenefit")
write.csv(cstat_res, "./res/counter_factual_cstat.csv", row.names = FALSE)
