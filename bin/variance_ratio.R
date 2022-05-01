# Script containin functions to perform Kurtosis adjusted variance ratio test for hetereogeneity
source("./bin/heterogeneity_presence_test.R")

trial_hte_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
                "NCT00460265", # CF takes a long time
                "NCT00041119_length", "NCT00041119_chemo",
                "NCT00003299", "NCT00119613")

var_res_df <- data.frame(trial_name = character(), stat = numeric(), pval = numeric())
for (trial in trial_hte_ls) {
    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
    
    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]
    # variance ratio test
    vrt_res <- variance_ratio_test(Y = Y, Z = W)
    var_res_df <- rbind(var_res_df, c(trial, vrt_res))
    
}
colnames(var_res_df) <- c("trial", "kurtosis_stat", "kurtosis_pval")
write.csv(var_res_df, "./res/variance_ratio_test.csv", row.names = FALSE)


## Omnibus test for systemic variations
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

    lm_df <- cbind(data.frame(tau = tau[,1]), X)
    taulm <- lm(tau ~ ., data = lm_df)
    stev <- try(wald.test(Sigma = vcov(taulm), b = coef(taulm), Terms = 2:dim(X)[2]))
    print(summary(taulm))
    if (stev != "try-error") {
        print(stev)
    } else {
        print("Can't perform Wald!")
    }

    # likelihood test by eliminating covaraites from model
    for (covar in colnames(X)) {
        #fit reduced model
        model_reduced <- lm(as.formula(paste0("tau ~ . - ", covar)), data = lm_df)

        #perform likelihood ratio test for differences in models
        print(lrtest(taulm, model_reduced))
    }

    # XGBoost to decompose effect modifier
    xgb_res <- cvboost(x = X, y = tau[,1], objective="reg:squarederror")
    shap_values <- shap.values(xgb_model = xgb_res, X_train = X)
}


# linearity check
pdf(paste("./res/plots/", trial, "_int_plot.pdf"))
interact_plot(taulm, pred = , modx = , plot.points = TRUE)
dev.off()

# Confidence interval band


# SHAP values
suppressPackageStartupMessages({
library("SHAPforxgboost"); library("ggplot2"); library("xgboost")
library("data.table"); library("here")
})
xgb_res <- readRDS("./dat/NCT00119613_xgb_obj.rds")
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
pdf(paste0("./res/plot/", trial, "_shap_plot.pdf"))
shap.plot.summary(shap_long)
dev.off()

# Interactions
# color feature: If "auto", will select the feature "c" minimizing the variance of the shap value given x and c, which can be viewed as a heuristic for the strongest interaction.
pdf(paste0("./res/plot/", trial, "_HGB_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'B_HGB', color_feature = 'auto') + ggtitle("(A) SHAP values of B_HGB vs. B_HGB")
dev.off()

pdf(paste0("./res/plot/", trial, "_AGE_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'AGE', color_feature = 'auto') + ggtitle("(A) SHAP values of AGE vs. AGE")
dev.off()
# not much of a pattern

# multiple interactions
pdf(paste0("./res/plot/", trial, "_multi-int_plot.pdf"))
fit_list <- lapply(names(shap_values$mean_shap_score)[1:6], 
                   shap.plot.dependence, data_long = shap_long)
gridExtra::grid.arrange(grobs = fit_list, ncol = 2)
dev.off()

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