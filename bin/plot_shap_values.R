source("./bin/load_lib.R")
# SHAP values
suppressPackageStartupMessages({
library("SHAPforxgboost"); library("ggplot2"); library("xgboost")
library("data.table"); library("here"); library("svglite"); library(broom)
})

min_mse_method <- read.csv("./res/best_tau_estimators.csv")
pos_report <- c("NCT00113763", "NCT00339183", "NCT00364013", "NCT00460265")
trial_choice <- filter(min_mse_method, mse < 150 & trial %in% pos_report)
colnames(trial_choice)[1] <- "trialID"


for (j in 1:dim(trial_choice)[1]) {

    trial <- trial_choice[j,1]
    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
    
    load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]

    outcome <- trial_choice[j,4]
    trial_best_method <- trial_choice[j,2]
    tau <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", outcome, "_", trial_best_method, "_tau_estimates.csv"))

    # Omnibus test for systemic variations
    cv_model <- cv.glmnet(X, tau[,1], keep = TRUE, family = "gaussian")
    best_lambda <- cv_model$lambda.min
    best_model <- glmnet(X, tau[,1], alpha = 1, lambda = best_lambda)
    lasso_coeff <- coef(best_model)
    lasso_coeff <- lasso_coeff[2:length(lasso_coeff)]
    X <- X[,lasso_coeff != 0]
    lm_df <- cbind(data.frame(tau = tau[,1]), X)
    taulm <- lm(tau ~ ., data = lm_df)
    # stev <- try(wald.test(Sigma = vcov(taulm), b = coef(taulm), Terms = 2:dim(X)[2]))
    print(summary(taulm))

    if (j == 1) {
        omnibus <- glance(taulm)
    } else {
        omnibus <- rbind(omnibus, glance(taulm))
    }

    # if (class(stev) != "try-error") {
    #     print(stev)
    # } else {
    #     print("Can't perform Wald!")
    # }

    sig_coeff <- c()
    # likelihood test by eliminating covaraites from model
    for (covar in colnames(X)) {
        #fit reduced model
        model_reduced <- lm(as.formula(paste0("tau ~ . - `", covar, "`")), data = lm_df)
        lrtres <- lrtest(taulm, model_reduced)
        #perform likelihood ratio test for differences in models
        if (covar == colnames(X)[1]) {
            lrt_df <- tidy(lrtres)
        } else {
            lrt_df <- rbind(lrt_df, tidy(lrtres)[2,])
        }
        if (lrtres$`Pr(>Chisq)`[2] < 0.05) {
            sig_coeff <- c(sig_coeff, covar)
        }
    }
    lrt_df <- cbind(c("full model", colnames(X)), lrt_df)
    colnames(lrt_df) <- c("Covariate", "DF", "LogLik", "diff", "statistic", "p_value")
    # write.csv(lrt_df, file = paste0("./res/omnibus/lrt/", trial, "_", outcome, "_", trial_best_method, "_LRT.csv"), row.names = FALSE)

    omnibus <- cbind(trial_choice[,c(1,4)], omnibus)
    # write.csv(omnibus, file = paste0("./res/omnibus/", "omnibus_res.csv"), row.names = FALSE)
    # XGBoost to decompose effect modifier
    xgb_res <- try(readRDS(paste0("./dat/xgb_model/", trial, "_", outcome, "_", trial_best_method, "_xgb_model.rds"))) # generated from cvboost(x = X, y = tau[,1], objective="reg:squarederror")
    if (class(xgb_res) == "try-error")
        next
    
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
    svglite(paste0("./res/shap_plots/", trial, "_", outcome, "_", trial_best_method, "_shap_plot.svg"))
    print(shap.plot.summary(shap_long))
    dev.off()
    # , height = ceiling(8 * length(sig_coeff) / 4)
    if (length(sig_coeff > 1)) {
        
        # multiple interactions
        svglite(paste0("./res/shap_plots/", trial, "_", outcome, "_", trial_best_method, "_multi-int_plot.svg"))
        fit_list <- lapply(sig_coeff, 
                        shap.plot.dependence, data_long = shap_long, color_feature = 'auto')
        gridExtra::grid.arrange(grobs = fit_list, ncol = 3, nrow = 3)
        dev.off()
    }
}

# Plot AUC
best_iter <- xgb_res$best_xgb_cvfit$evaluation_log

library(pROC)
plot(pROC::roc(response = tau[,1],
               predictor = predict(mod, newdata = X),
               levels=c(0, 1)), lwd=1.5)

explainer = buildExplainer(xgb.model,xgb.train.data, type="binary", base_score = 0.5, trees_idx = NULL)
pred.breakdown = explainPredictions(xgb.model, explainer, xgb.test.data)




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
svglite(paste0("./res/plot/", trial, "_KRASe2_int_plot.svg"))
shap.plot.dependence(data_long = shap_long, x = 'KRAS_exon_2_(c12/13)', color_feature = 'auto') + ggtitle("SHAP values of KRAS_exon_2_(c12/13) vs. KRAS_exon_2_(c12/13)")
dev.off()

svglite(paste0("./res/plot/", trial, "_BWEIGHT_int_plot.svg"))
shap.plot.dependence(data_long = shap_long, x = 'B_WEIGHT', color_feature = 'auto') + ggtitle("SHAP values of Baseline weight vs. Baseline weight Performance Status")
dev.off()

pdf(paste0("./res/plot/", trial, "_AGE_int_plot.pdf"))
shap.plot.dependence(data_long = shap_long, x = 'AGE', color_feature = 'auto') + ggtitle("SHAP values of age vs. age")
dev.off()

svglite(paste0("./res/plot/", trial, "_KRASe4_int_plot.svg"))
shap.plot.dependence(data_long = shap_long, x = 'KRAS_exon_4_(c117/146)', color_feature = 'auto') + ggtitle("SHAP values of KRAS_exon_4_(c117/146) vs. KRAS_exon_4_(c117/146)")
dev.off()

