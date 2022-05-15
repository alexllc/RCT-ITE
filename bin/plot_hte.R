#' Script to plot tau value distributions
source("./bin/load_lib.R")
library(ggpubr)

# trial_list <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
#                 "NCT00460265", # CF takes a long time
#                 "NCT00041119_length", "NCT00041119_chemo",
#                 "NCT00003299", "NCT00119613")

trial_list <- c("NCT00460265", "NCT00364013", "NCT00113763")

min_mse_method <- read.csv("./res/crossfit_rloss/best_tau_estimators.csv")

for (trial in trial_list) {

    tau_method <- min_mse_method[min_mse_method$trial == trial, "best_tau_method"]
    tau_values <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", tau_method , "_tau_estimates.csv"))

    if (tau_method == "CF") {
        tau <- tau_values$predictions
    } else {
        tau <- tau_values[,1]
    }

    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]
    # plot distribution of tau values

    plot_df <- data.frame(SUBJID = 1:length(tau), tau = tau, sign = ifelse(W, "Treatment", "Control"))

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
            legend.title = "TE"
            ) + theme(axis.text.x=element_blank()) + geom_hline(yintercept=1.6,         color = "red", size=2)
            )
    dev.off()
}