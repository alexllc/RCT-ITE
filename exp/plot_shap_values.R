source("./bin/load_lib.R")
# SHAP values
suppressPackageStartupMessages({
library("SHAPforxgboost"); library("ggplot2"); library("xgboost")
library("data.table"); library("here"); library("svglite"); library(broom)
})

ite_list <- list.files("./res/ite_tau_estimates/")
ite_list <- ite_list[grep("NCT*", ite_list)]

for (ite_model in ite_list) {

    trial <- strsplit(ite_model, "_")[[1]][1]
    outcome <- strsplit(ite_model, "_")[[1]][2]
    trial_best_method <- strsplit(ite_model, "_")[[1]][3]

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

    tau <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", outcome, "_", trial_best_method, "_tau_estimates.csv"))

    # Omnibus test for systemic variations
    cv_model <- cv.glmnet(X, tau[,1], keep = TRUE, family = "gaussian")
    best_lambda <- cv_model$lambda.min
    best_model <- glmnet(X, tau[,1], alpha = 1, lambda = best_lambda)
    lasso_coeff <- coef(best_model)
    lasso_coeff <- lasso_coeff[2:length(lasso_coeff)]
    lm_X <- X[,lasso_coeff != 0]
    lm_df <- cbind(data.frame(tau = tau[,1]), lm_X)
    taulm <- lm(tau ~ ., data = lm_df)
    stev <- try(wald.test(Sigma = vcov(taulm), b = coef(taulm), Terms = 2:dim(lm_X)[2]))
    print(summary(taulm))

    if (class(stev) != "try-error") {
        print(stev)
        omni_res <- cbind(trial, outcome, glance(taulm), t(as.data.frame(stev$result$chi2)))
    } else {
        print("Can't perform Wald!")
        stev <- rep(NA, 3)
        names(stev) <- c("chi2", "df", "P")
        omni_res <- cbind(trial, outcome, glance(taulm), t(as.data.frame(stev)))

    }
    if (ite_model == "NCT00052910_OS_CBoost_tau_estimates.csv") {
        omnibus <- omni_res
    } else {
        omnibus <- rbind(omnibus, omni_res)
    }

    # write.csv(omnibus, file = "./res/omnibus/omnibus_lr_res.csv", row.names = FALSE)

    sig_coeff <- c()
    # likelihood test by eliminating covaraites from model
    for (covar in colnames(lm_X)) {
        #fit reduced model
        model_reduced <- lm(as.formula(paste0("tau ~ . - `", covar, "`")), data = lm_df)
        lrtres <- lrtest(taulm, model_reduced)
        #perform likelihood ratio test for differences in models
        if (covar == colnames(lm_X)[1]) {
            lrt_df <- tidy(lrtres)
        } else {
            lrt_df <- rbind(lrt_df, tidy(lrtres)[2,])
        }
        if (lrtres$`Pr(>Chisq)`[2] < 0.05) {
            sig_coeff <- c(sig_coeff, covar)
        }
    }
    lrt_df <- cbind(c("full model", colnames(lm_X)), lrt_df)
    colnames(lrt_df) <- c("Covariate", "DF", "LogLik", "diff", "statistic", "p_value")
    write.csv(lrt_df, file = paste0("./res/omnibus/lrt/", trial, "_", outcome, "_", trial_best_method, "_LRT.csv"), row.names = FALSE)

    omnibus <- cbind(trial_choice[,c(1,4)], omnibus)
    # write.csv(omnibus, file = paste0("./res/omnibus/", "omnibus_res.csv"), row.names = FALSE)
    # XGBoost to decompose effect modifier
    xgb_res <- try(readRDS(paste0("./dat/xgb_model/", trial, "_", outcome, "_", trial_best_method, "_xgb_model.rds"))) # generated from cvboost(x = X, y = tau[,1], objective="reg:squarederror")
    if (class(xgb_res) == "try-error")
        next
    
    mod <- xgb_res$xgb_fit
    X <- as.matrix(get(trial)[[1]])
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
    
    plot_coeff <- filter(lrt_df, p_value < 0.05)
    plot_coeff <- plot_coeff$Covariate

    for (clinvar in plot_coeff) {
        if (outcome != "ORR" & outcome != "RSP") {
            shap_thresh <- filter(shap_long, variable == clinvar & value > 0)
            dir_lm <- lm(shap_thresh$value ~ shap_thresh$rfvalue)
            dir_lm <- tidy(dir_lm)
            if (dir_lm$estimate[2] > 0) {
                thresh <- min(shap_thresh$rfvalue)
                print(clinvar)
                print(paste0("Min: ", thresh))
            } else {
                thresh <- max(shap_thresh$rfvalue)
                print(clinvar)
                print(paste0("Max: ", thresh))
            }
        } else {
            shap_thresh <- filter(shap_long, variable == clinvar & value < 0)
            dir_lm <- lm(shap_thresh$value ~ shap_thresh$rfvalue)
            dir_lm <- tidy(dir_lm)
            if (dir_lm$estimate[2] < 0) {
                thresh <- min(shap_thresh$rfvalue)
                print(clinvar)
                print(paste0("Min: ", thresh))
            } else {
                thresh <- max(shap_thresh$rfvalue)
                print(clinvar)
                print(paste0("Max: ", thresh))
            }
        }
    }

    # multiple interactions
    for (k in 1:ceiling(length(plot_coeff) / 4)) {
        svglite(paste0("./res/shap_plots/", trial, "_", outcome, "_", trial_best_method, "_multi-int_plot_", k,".svg"))
        if (k == 1)
            start_counter <- 1
        selection <- plot_coeff[start_counter:ifelse((start_counter+3) > length(plot_coeff),  length(plot_coeff), (start_counter+3))] #
        fit_list <- lapply(selection, 
                        shap.plot.dependence, data_long = shap_long, color_feature = 'auto', size0 = 1)
        gridExtra::grid.arrange(grobs = fit_list, ncol = 2, nrow = 2)
        dev.off()
        start_counter <- start_counter + 4
        print(selection)
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

