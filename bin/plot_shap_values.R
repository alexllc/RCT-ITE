source("./bin/load_lib.R")
# SHAP values
suppressPackageStartupMessages({
library("SHAPforxgboost"); library("ggplot2"); library("xgboost")
library("data.table"); library("here")
})

trial_hte_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",  "NCT00460265", "NCT00041119_length", "NCT00041119_chemo", "NCT00003299", "NCT00119613")

min_mse_method <- read.csv("./res/crossfit_rloss/best_tau_estimators.csv")


for (trial in trial_hte_ls) {
    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
    
    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]

    tau <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", min_mse_method[min_mse_method$trial == trial, "best_tau_method"], "_tau_estimates.csv"))

    # XGBoost to decompose effect modifier
    xgb_res <- readRDS(paste0("./dat/xgb_model/", trial, "_xgb_model.rds")) # generated from cvboost(x = X, y = tau[,1], objective="reg:squarederror")
    
    mod <- xgb_res$xgb_fit
    shap_values <- shap.values(xgb_model = mod, X_train = X)
    shap_values$mean_shap_score

    shap_data <- copy(shap_values$shap_score)
    shap_data[, BIAS := shap_values$BIAS0]
    pred_mod <- predict(mod, X, iteration_range = 10)
    shap_data[, `:=`(rowSum = round(rowSums(shap_data),6), pred_mod = round(pred_mod,6))]

    # To prepare the long-format data:
    shap_long <- shap.prep(xgb_model = mod, X_train = X)
    # is the same as: using given shap_contrib
    shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = X)

    # summary plot
    pdf(paste0("./res/plot/", trial, "/", trial, "_shap_plot.pdf"))
    print(shap.plot.summary(shap_long))
    dev.off()

    png(paste0("./res/plot/", trial, "/", trial, "_shap_plot.png"))
    print(shap.plot.summary(shap_long))
    dev.off()

    # multiple interactions
    pdf(paste0("./res/plot/", trial, "/", trial, "_multi-int_plot.pdf"))
    fit_list <- lapply(names(shap_values$mean_shap_score)[1:6], 
                    shap.plot.dependence, data_long = shap_long)
    gridExtra::grid.arrange(grobs = fit_list, ncol = 2)
    dev.off()

    png(paste0("./res/plot/", trial, "/", trial, "_multi-int_plot.png"))
    fit_list <- lapply(names(shap_values$mean_shap_score)[1:6], 
                    shap.plot.dependence, data_long = shap_long)
    gridExtra::grid.arrange(grobs = fit_list, ncol = 2)
    dev.off()
}

#############################################################################
## NCT00119613
#############################################################################

# Interactions
# color feature: If "auto", will select the feature "c" minimizing the variance of the shap value given x and c, which can be viewed as a heuristic for the strongest interaction.
pdf(paste0("./res/plot/", trial, "_HGB_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'B_HGB', color_feature = 'auto') + ggtitle("(A) SHAP values of B_HGB vs. B_HGB")
dev.off()

pdf(paste0("./res/plot/", trial, "_AGE_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'AGE', color_feature = 'auto') + ggtitle("(A) SHAP values of AGE vs. AGE")
dev.off()
# not much of a pattern


# Separating interaction effects from main effects
shap_int <- shap.prep.interaction(xgb_mod = mod, X_train = X)

g3 <- shap.plot.dependence(data_long = shap_long,
                           data_int = shap_int,
                           x= "AGE", y = "AGE", 
                           color_feature = "auto")
g4 <- shap.plot.dependence(data_long = shap_long,
                           data_int = shap_int,
                           x= "B_HGB", y = "B_HGB", 
                           color_feature = "auto")

pdf(paste0("./res/plot/", trial, "_multi-int_plot.pdf"))          
gridExtra::grid.arrange(g3, g4, ncol=2)
dev.off()

#############################################################################
# NCT00113763
#############################################################################

# Interactions
pdf(paste0("./res/plot/", trial, "_KRASe2_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'KRAS_exon_2_(c12/13)', color_feature = 'auto') + ggtitle("SHAP values of KRAS_exon_2_(c12/13) vs. KRAS_exon_2_(c12/13)")
dev.off()

pdf(paste0("./res/plot/", trial, "_BECOG_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'B_ECOG', color_feature = 'auto') + ggtitle("SHAP values of Baseline ECOG Performance Status
 vs. Baseline ECOG Performance Status")
dev.off()

pdf(paste0("./res/plot/", trial, "_AGE_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'AGE', color_feature = 'auto') + ggtitle("SHAP values of age vs. age")
dev.off()

pdf(paste0("./res/plot/", trial, "_KRASe4_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'KRAS_exon_4_(c117/146)', color_feature = 'auto') + ggtitle("SHAP values of KRAS_exon_4_(c117/146) vs. KRAS_exon_4_(c117/146)")
dev.off()

