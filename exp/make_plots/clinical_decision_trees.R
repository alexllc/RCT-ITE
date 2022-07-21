source("./bin/load_lib.R")
source("./exp/make_plots/plot_tree_func.r")
suppressPackageStartupMessages({
library("SHAPforxgboost"); library("ggplot2"); library("xgboost")
library("data.table"); library("here"); 
library(xgboostExplainer); library(rpart.plot) ; library(rattle); library(partykit)
})


sig_hte_trials <- c("NCT00113763")

# Find sig. HTE trial and outcome pairs

calib_sum_files <- list.files("./res/calib/backup/calib_avg/")
calib_sum_files <- calib_sum_files[grep("*.csv", calib_sum_files)]

sig_trial_outcome <- data.frame()

for (sum_file in calib_sum_files) {
    calib_res <- read.csv(paste0("./res/calib/backup/calib_avg/", sum_file))
    sig_res <- filter(calib_res, tau_method == "Rstack" & hte.pval < 0.1)
    if (dim(sig_res)[1] == 0) {
        next
    } else {
        trial <- unlist(lapply(strsplit(sum_file, "_|\\."), function(X) X[1]))
        sig_trial_outcome <- rbind(sig_trial_outcome, cbind(trial, sig_res$outcome))
    }
}


for (i in 1:dim(sig_trial_outcome)[1]) {

    trial <- sig_trial_outcome[i,1]
    outcome <- sig_trial_outcome[i,2]

    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }

    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]

    Y_list <- get(paste0(outcome, "_Y_list"))

    imp_type_Y <- "efronYn"
    Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

    if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
        Y <- as.numeric(Y_list[[1]])
        imp_type_Y <- "efron"
    }
    rmv_id <- check_lm_rmv(X = X, Y = Y)

    X <- X[-rmv_id,]
    Y <- Y[-rmv_id]
    W <- W[-rmv_id]

    tau <- try(read.csv(paste0("./res/ite_tau_estimates/", trial, "_", imp_type_Y, "_", outcome, "_Rstack_tau_estimates.csv")))

    old_tau <- read.csv("./res_PFS_only/ite_tau_estimates/NCT00041119_length_CMARS_tau_estimates.csv")

    if (class(tau) == "try-error") {
        next
    }

    tau_vec <- old_tau[,1]

    tree_df <- data.frame(tau_vec, X)

    # rf <- randomForest(tau_vec ~., data = tree_df, ntree=5000, maxnodes = 10)
    # best_tree <- reprtree::ReprTree(rf, newdata = X)
    # best_tree_id <- as.numeric(names(best_tree$trees))
    # best_tree_extract <- getTree(rf, k = best_tree_id)

    # pdf("best_tree.pdf")
    # tree_func(rf, best_tree_id)
    # dev.off()

    rtree_model <- rpart(tau_vec~., data=tree_df, control=rpart.control(maxdepth=4))

    pdf(paste0("./clin_decision_tree/", trial, "_", outcome, "_tree_plot.pdf"))
    fancyRpartPlot(rtree_model) 
    dev.off()

}