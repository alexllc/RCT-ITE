#' Script to plot tau value distributions
source("./bin/load_lib.R")
library(ggpubr)

# trial_list <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
#                 "NCT00460265", # CF takes a long time
#                 "NCT00041119_length", "NCT00041119_chemo",
#                 "NCT00003299", "NCT00119613")

# SIG TRIALS only
#  trial_list <- c("NCT00460265", "NCT00364013", "NCT00113763")

min_mse_method <- read.csv("./res/best_tau_estimators.csv")
trial_choice <- filter(min_mse_method, mse < 150)
colnames(trial_choice)[1] <- "trialID"

tau_summary <- data.frame(trial = character(), Min = numeric(), FirstQu = numeric(), Median  = numeric(), Mean  = numeric(), ThirdQu = numeric(), Max = numeric())

for (trial in trial_list) {

    tau_method <- min_mse_method[min_mse_method$trial == trial, "best_tau_method"]
    tau_values <- read.csv(paste0("./res_PFS_only/ite_tau_estimates/", trial, "_", tau_method , "_tau_estimates.csv"))

    if (tau_method == "CF") {
        tau <- tau_values$predictions
    } else {
        tau <- tau_values[,1]
    }
    tau_summary <- rbind(tau_summary,  c(summary(tau), var(tau)))
}
tau_summary <- cbind(trial_list, tau_summary)

colnames(tau_summary) <-  c("Trial", "Min", "1stQu", "Median", "Mean", "3rdQu", "Max", "Var")

    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]
    # plot distribution of tau values

    plot_df <- data.frame(SUBJID = 1:length(tau), tau = tau, sign = ifelse(W, "Intervention", "Control"))

    pdf(file = paste0(trial, "_", tau_method, "_plot.pdf"), width = 30, height = 20)
    print(ggbarplot(plot_df, x = "SUBJID", y = "tau",
            fill = "sign",           # change fill color by mpg_level
            color = "white",            # Set bar border colors to white
            palette = "jco",            # jco journal color palett. see ?ggpar
            sort.val = "asc",           # Sort the value in ascending order
            sort.by.groups = TRUE,     # Don't sort inside each group
            x.text.angle = 90,          # Rotate vertically x axis texts
            ylab = trial,
            xlab = FALSE,
            legend.title = "Treatment assignment"
            ) + theme(axis.text.x=element_blank()) + geom_hline(yintercept=mean(tau),         color = "red", size=2)
            )
    dev.off()
}

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