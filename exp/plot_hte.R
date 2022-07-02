#' Script to plot tau value distributions
source("./bin/load_lib.R")
library(ggpubr)

# trial_list <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
#                 "NCT00460265", # CF takes a long time
#                 "NCT00041119_length", "NCT00041119_chemo",
#                 "NCT00003299", "NCT00119613")

# SIG TRIALS only
#  trial_list <- c("NCT00460265", "NCT00364013", "NCT00113763")

#
# Save summary statistics of tau estimates
#
min_mse_method <- read.csv("./res/best_tau_estimators.csv")
pos_report <- c("NCT00113763", "NCT00115765", "NCT00339183", "NCT00364013", "NCT00460265")
trial_choice <- filter(min_mse_method, trial %in% pos_report)
colnames(trial_choice)[1] <- "trialID"

tau_summary <- data.frame(trial = character(), Min = numeric(), FirstQu = numeric(), Median  = numeric(), Mean  = numeric(), ThirdQu = numeric(), Max = numeric())

for (i in 1:dim(trial_choice)[1]) {
    trial <- trial_choice[i,1]
    tau_method <- trial_choice[i,2]
    outcome <- trial_choice[i,4]
    tau_values <- try(read.csv(paste0("./res/ite_tau_estimates/", trial, "_", outcome, "_", tau_method , "_tau_estimates.csv")))
    if (class(tau_values) == "try-error") {
        next
    }
    tau <- tau_values[,1]
    tau_summary <- rbind(tau_summary, c(trial, outcome, summary(tau), var(tau)))

    colnames(tau_summary) <-  c("Trial",  "outcome","Min", "1stQu", "Median", "Mean", "3rdQu", "Max", "Var")
    for(k in 1:dim(tau_summary)[2]) {
        if (!(k == 1 | k == 2))
            tau_summary[,k] <- as.numeric(tau_summary[,k])
    }
    write.csv(tau_summary, "./res/ite_tau_estimates/ITE_summary.csv", row.names = FALSE)
}

trial_units <- read.csv("./res/units.csv")
colnames(trial_units)[2] <- "outcome_type" #1:dim(trial_choice)[1]
for (i in 7:8) {
    trial <- trial_choice[i,1]
    tau_method <- trial_choice[i,2]
    outcome <- trial_choice[i,4]
    tau_values <- try(read.csv(paste0("./res/ite_tau_estimates/", trial, "_", outcome, "_", tau_method , "_tau_estimates.csv")))
    if (class(tau_values) == "try-error") {
        next
    }
    tau <- tau_values[,1]

    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]


    message(paste0("Processing outcome: ", outcome))

    Y_list <- get(paste0(outcome, "_Y_list"))

    imp_type_Y <- "efronYn"
    Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

    if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
        Y <- as.numeric(Y_list[[1]])
        imp_type_Y <- "efron"
    }
    
    # plot distribution of tau values
    outcome_unit <- filter(trial_units, outcome_type == outcome, Trial == trial)
    outcome_unit <- outcome_unit$unit
    
    plot_df <- data.frame(SUBJID = 1:length(tau), tau = tau, sign = ifelse(W, "Intervention", "Control"))
    png(file = paste0("./res/ite_tau_estimates/plots/", trial, "_", outcome, "_", tau_method, "_KRAS_stratify_plot.png"), width = 2560, height = 1600)
    print(ggbarplot(plot_df, x = "SUBJID", y = "tau",
            fill = "sign",           # change fill color by mpg_level
            color = "white",            # Set bar border colors to white
            palette = "jco",            # jco journal color palett. see ?ggpar
            sort.val = "asc",           # Sort the value in ascending order
            sort.by.groups = TRUE,     # Don't sort inside each group
            x.text.angle = 90,          # Rotate vertically x axis texts
            xaxt = "n",
            ylab = outcome_unit,
            xlab = FALSE,
            legend.title = "Original treatment assignment"
            ) + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 40))
              + geom_hline(yintercept=mean(tau), color = "red", size=2)
            )
    dev.off()
}

#
# Save CB parameters used
#

first_cb <- cb_param

for (trial in trial_list) {
     cb_param_tune_res <- read.csv(paste0("./dat/cb_param_PFS_only/single_outcome/", trial, "_cb_param.csv"))
     if ("X" %in% colnames(cb_param_tune_res)) {
         cb_param_tune_res$X <- NULL
     }
     if ("nfolds" %in% colnames(cb_param_tune_res)) {
         cb_param_tune_res$nfolds <- NULL
     }
     cb_param_tune_res <- filter(cb_param_tune_res, num_trees_min_effect != 1) # avoid using parameteres with only 1 boosting tree
    cb_param <- cb_param_tune_res[which.min(cb_param_tune_res$mean_cvm_effect),]
    first_cb <- rbind(first_cb, cb_param)
}
first_cb <- cbind(trial_list, first_cb)
write.csv(first_cb, file = "./dat/cb_param_PFS_only/cb_param_chosen.csv", row.names = FALSE)


#
# Treatment discordance ratio
#
# When the estimated counterfactual treatment effect is negative but the 