suppressPackageStartupMessages({
library(tidyr);
library(data.table); # you should load tidyr before loading dplyr
library(dplyr); library(impute); library(naniar);
# Survival analysis
library(survival); library(survminer); library(imputeYn); library(pseudo)

# Causal inference
library(causalToolbox) # X-learner
library(grf); library(causalLearning) # Powers

# R-learner
library(rlearner); library(KRLS2); library(glmnet); library(nnls); library(caret) # cross validation 
# Find HTE
library(hettx); library(formula.tools) ; # a bug in the package failed to import forumla.tools s.t. FRTCI cannot get variables from the formula

# post hoc analysis
library(aod); library(lmtest);

library(tidyverse); library(Hmisc); library(corrplot); # for removing colinearity

library(broom); library(ggfortify)
})

source("./bin/impute_survival.R")

na_strings <- c("", "NA", "98", 98, "99", 99, "missing", "Unknown", "unknown", "<NA>", "Unable to evaluate", "Non Evaluable", "Failure", "failure", "Not applicable", "Not Applicable", "Not done", "Not Done","Other", "other")

coalesce_by_column <- function(df) {
  return(coalesce(df[1], df[2]))
}


nnls <- function(M, v, constrained) {
    Dmat <- t(M) %*% M
    dvec <- t(M) %*% v
    Amat <- matrix(0, ncol(M), sum(constrained))
    cons <- which(constrained)
    for(iter in 1:length(cons)) {
        Amat[cons[iter], iter] <- 1
    }
    bvec <- rep(0, sum(constrained))
    soln <- quadprog::solve.QP(Dmat, dvec, Amat, bvec) 
    soln$solution
}

rloss <- function(tau.pred = NULL, Y = NULL, W = NULL, Y.hat = NULL, prob = NULL) {
    RESP = Y - Y.hat
    R.mat = cbind(1, W - prob,
                    (W - prob) * tau.pred)

    learn_coeff = nnls(R.mat, RESP, constrained = c(FALSE, FALSE, TRUE))

    print("coefs")
    print(learn_coeff)
    debiased_tau <- learn_coeff[2] + learn_coeff[3] * tau.pred
    mse_debiased <- mean((Y - Y.hat - debiased_tau * (W - prob))^2)
    mse <- mean((Y - Y.hat - tau.pred * (W - prob))^2)

    return(list(nnls_coeff = learn_coeff, debiased_tau = debiased_tau, mse_debiased = mse_debiased, mse = mse))
}

#
# data cleaning functions
#

cat2num <- function(df, verbose = TRUE) {
    col_class <- unlist(lapply(df, function(X) class(X)))
    cat_col <- which(col_class == "character")
    ddt <- data.frame(cat_variable = cbind(names(cat_col)), cat_levels = seq_along(cat_col))
    i <- 1
    for (cat_indx in cat_col) {
        factor_col <- as.factor(df[,cat_indx])
        save_lv <- sapply(levels(factor_col), function(X) gsub("$", "=", X))
        df[,cat_indx] <- as.numeric(factor_col)
        ddt[i, 2] <- paste0(apply(rbind(save_lv, seq_along(save_lv)), 2, function(X) paste0(X, collapse = "")), collapse = "|")
        i <- i+1
    }
    return(list(df, ddt))
}


scaleNimpute <- function(df) {
    df <- scale(df, center = TRUE, scale = TRUE)
    df <- knn
    df <- t(apply(df, 1, function(r) r * attr(df,'scaled:scale') + attr(df, 'scaled:center')))
}


hdf <- function(df, n = 5){
    print(head(as.data.frame(df), n = n))
}


impute_df_missing  <- function(clin_df = NULL, save_ddt = FALSE) {

    # convert categorical to numeric
    num_df <- cat2num(clin_df)
    num_clin_df <- num_df[[1]]
    clin_df_ddt <- num_df[[2]]

    #
    # impute missing data if necessary
    #
    if (any(is.na(num_clin_df))) {
        pre_proc <- scale(num_clin_df, center = TRUE, scale = TRUE)
        imp_clin_df <- impute.knn(pre_proc, k=10)
        message(c("Seed used: ", imp_clin_df$rng.seed))
        for (i in 1:dim(imp_clin_df$data)[2]) {
            imp_clin_df$data[,i] <- imp_clin_df$data[,i] * attr(pre_proc, "scaled:scale")[i] +  attr(pre_proc, "scaled:center")[i]
        }
    }

    # check if KNN imputation is successful
    message("Is there still missing data? ", any(is.na(imp_clin_df$data)))
    imp_clin_df <- as.data.frame(imp_clin_df$data)

    if(save_ddt)
        return(list(imp_clin_df, num_clin_df))
    else
        return(imp_clin_df)
}

missing_too_much <- function(clin_df, missing_threshold = 0.3) {

    for(col in colnames(clin_df)) {
        missing_prop <- sum(is.na(clin_df[[col]])) / dim(clin_df)[1]
        zero_prop <- sum(clin_df[[col]] == 0, na.rm = TRUE) / dim(clin_df)[1]
        if (missing_prop > missing_threshold) {
            print(paste0("More than 30% missing data! Removing from df.: ", col))
            clin_df[[col]] <- NULL
        } else if(length(unique(na.omit(clin_df[[col]]))) == 1 | var(na.omit(clin_df[[col]]) == 0)) {
                print(paste0("All values are the same, removing: ", col))
                clin_df[[col]] <- NULL
        } else if (zero_prop > 0.8) {
            print(paste0("More than 80% of values are 0s, removing: ", col))
            clin_df[[col]] <- NULL
        }
    }

    return(clin_df)
}

hdf <- function(df, n = 10) {
    print(head(as.data.frame(df), n = n))
}


tau_calibration <- function(tau.hat = NULL, Y = NULL, m.hat = NULL, W = NULL, W.hat = NULL, vcov.type = "HC3") {
    ate <- mean(Y[which(W == 1)]) - mean(Y[which(W == 0)])
    DF <- data.frame(target = Y - m.hat, mean_pred = (W - W.hat) * ate, hetereogeneity = (W - W.hat) * tau.hat)

    best.linear.predictor <-
    lm(target ~ mean_pred + hetereogeneity + 0,
      data = DF
    )
    blp.summary <- lmtest::coeftest(best.linear.predictor,
        vcov = sandwich::vcovHC(best.linear.predictor, type = vcov.type)
    )
    attr(blp.summary, "method") <-
        paste0("Best linear fit using ITE estimations (on held-out data)\n",
        "as well as the ATE as regressors, along\n",
        "with one-sided heteroskedasticity-robust ",
        "(", vcov.type, ") SEs"
        )
    # convert to one-sided p-values
    dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
    blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[, 4] / 2, blp.summary[, 4] / 2)
    print(blp.summary)

    return(tidy(blp.summary))
}

# Function to determine whether to use boosting or LASSO to derive mhat
get_mhat <- function(X_train, W_train, Y_train, X_ho, Y_ho, binary_Y = FALSE, newdat = TRUE, fast = FALSE) {

    if (fast) {
        if (binary_Y) {
            Y_lasso <- cv.glmnet(X_train, Y_train, keep = TRUE, family = "binomial")
        } else {
            Y_lasso <- cv.glmnet(X_train, Y_train, keep = TRUE, family = "gaussian")
        }
        Y_hat_lasso <- predict(Y_lasso, newx = X_ho)[,1]
        print(RMSE(Y_hat_lasso, Y_ho))
        return(Y_hat_lasso)
    } else {
        Y_boost <- cvboost(X_train, Y_train, objective = "reg:squarederror", nthread = 40) # non-binary outcome
        Y_hat_boost <- predict(Y_boost, newx = X_ho)
        
        if (binary_Y) {
            Y_lasso <- cv.glmnet(X_train, Y_train, keep = TRUE, family = "binomial")
        } else {
            Y_lasso <- cv.glmnet(X_train, Y_train, keep = TRUE, family = "gaussian")
        }
        Y_hat_lasso <- predict(Y_lasso, newx = X_ho)[,1]

        # which method has a smaller CV error?
        print("RMSE of m hat estimated by boosting vs LASSO:")
        print(round(c(RMSE(Y_hat_boost, Y_ho), RMSE(Y_hat_lasso, Y_ho)), 4))
        
        if (!newdat) {
            Y_hat_boost <- predict(Y_boost)
            Y_hat_lasso <- predict(Y_lasso, newx = X_train)[,1]
            if (RMSE(Y_hat_boost, Y_ho) < RMSE(Y_hat_lasso, Y_ho)) {
                return(Y_hat_boost)
            } else {
                return(Y_hat_lasso)
            }
        } else {
            if (RMSE(Y_hat_boost, Y_ho) < RMSE(Y_hat_lasso, Y_ho)) {
                return(Y_hat_boost)
            } else {
                return(Y_hat_lasso)
            }
        }
    }
}

check_lm_rmv <- function(X = NULL, Y = NULL) {
    check_df <- data.frame(outcome = Y, trainX = X)
    model <- lm(outcome ~ ., data = check_df)
    print(vif_res <- car::vif(model))

    model.diag.metrics <- augment(model)

    model.diag.metrics <- model.diag.metrics %>% mutate(index = 1:nrow(model.diag.metrics)) 
  
    # Inspect the data
    head(model.diag.metrics)

    # outliers
    outlier_id <- filter(model.diag.metrics, abs(.std.resid) >= 3)
    outlier_id <- outlier_id$index

    # high leverage
    high_lev_thresh <- 2*(dim(X)[2] + 1)/ dim(X)[1]
    highlev_id <- filter(model.diag.metrics, .hat > high_lev_thresh)
    highlev_id <- highlev_id$index


    # high influence
    high_influ_thresh <- 4/(dim(X)[1] - dim(X)[2] - 1)
    highinflu_id <- filter(model.diag.metrics, .cooksd > high_influ_thresh)
    highinflu_id <- highinflu_id$index

    rmv_id <- unique(c(outlier_id, highlev_id, highinflu_id))
    return(rmv_id)
}