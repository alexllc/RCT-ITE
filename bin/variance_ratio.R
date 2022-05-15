# Script containin functions to perform Kurtosis adjusted variance ratio test for hetereogeneity
source("./bin/load_lib.R")
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

    sig_coeff <- c()
    # likelihood test by eliminating covaraites from model
    for (covar in colnames(X)) {
        #fit reduced model
        model_reduced <- lm(as.formula(paste0("tau ~ . - `", covar, "`")), data = lm_df)
        lrtres <- lrtest(taulm, model_reduced)
        #perform likelihood ratio test for differences in models
        print(lrtres)
        if (lrtres$`Pr(>Chisq)`[2] < 0.05) {
            sig_coeff <- c(sig_coeff, covar)
        }
    }


}
