# Function to find the best tau predictor using R-loss criteria averaging over several folds
find_best_tau_estimator <- funciton(Y = NULL, X = NULL, W = NULL, prob = 0.5, Q = 4, tuned_cb_param = FALSE) {

    # For testing only:
    source("./bin/load_CRC_RCT/load_NCT00339183.R")
    X <- as.matrix(NCT00339183[[1]])
    Y <- NCT00339183[[2]][[1]]
    W <- NCT00339183[[3]]

    # set up empty dataframe
    compare_rloss <- data.frame(b = numeric(), c = numeric(), alpha = numeric(), mse = numeric(), debiased_mse = numeric())

    # We need to take out 1/Q sample as the holdout test data, the rest of the samples are used for training
    holdout_sample <- sample(1:dim(X)[1], dim(X)[1] / Q)
    train_sample <- seq(1:dim(X)[1])[!(seq(1:dim(X)[1]) %in% holdout_sample)]
    ho_X <- X[holdout_sample,]
    ho_W <- W[holdout_sample]
    ho_Y <- Y[holdout_sample]

    train_X <- X[train_sample,]
    train_W <- W[train_sample]
    train_Y <- Y[train_sample]

    # Since CausalBoost is not implemented w.r.t. the causal forest in Imbens and Athey, there is not debiased.error to take out of context

    # let's first fit cvlasso or cv boosting to estimate nuisance component marginal response model (m.hat)

    Y.boost <- cvboost(train_X, train_Y, objective = "reg:squarederror", nthread = 40) # non-binary outcome
    Y.hat.boost <- predict(Y.boost, newx = ho_X)

    Y.lasso <- cv.glmnet(train_X, train_Y, keep = TRUE, family = "gaussian")
    Y.hat.lasso <- predict(Y.lasso, newx = ho_X)[,1]

    # which method has a smaller CV error?
    print(round(c(RMSE(Y.hat.boost, ho_Y), RMSE(Y.hat.lasso, ho_Y)), 4))
    
    if( RMSE(Y.hat.boost, ho_Y) < RMSE(Y.hat.lasso, ho_Y) ) {
        Y.hat <- Y.hat.boost
    } else {
        Y.hat <- Y.hat.lasso
    }
    #
    # Causal forest
    #

    # R-loss already implemented by the package and is contained in the debiased.error column in predict.causal_forest(). You can obtain MSE by mean(cf$debiased.error), the values were already squared.
    # since causal trees are honest already, there's no need to further split into k-fold estimation

    cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
    cf_tau_pred <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)

    cf_mse <- mean((ho_Y - Y.hat - cf_tau_pred$predictions * (ho_W - prob))^2)
    cf_rloss_res <- rloss(tau.pred = cf_tau_pred$predictions, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)

    compare_rloss <- rbind(compare_rloss, c(cf_rloss_res$nnls_coeff, cf_rloss_res$mse, cf_rloss_res$mse_debiased))


    ##################################################
    # Powers: CausalBoost, CausalMARS, PTOForest
    ##################################################

    if(tuned_cb_param) {
        param_list <- read.csv("./dat/NCT00339183_cb_param.csv")
        param_list <- param_list[which(param_list$mean_cvm_effect == min(param_list$mean_cvm_effect)),]
        param_list <- list(num.trees = param_list$num_trees_min_effect, maxleaves = param_list$max_leaves, eps = param_list$eps, splitSpread = param_list$split_spread)
    } else { # use default setting
        param_list <- list(num.trees = 500, maxleaves = 4, eps = 0.01, splitSpread = 0.1)
    }
    # now we estimate tau from CausalBoosting
    
    cb <- do.call(causalBoosting, append(list(x = train_X, tx = train_W, y = train_Y), as.list(param_list)))
    tau_cb_pred <- predict(cb, newx = ho_X)
    tau_cb_pred <- tau_cb_pred[,param_list$num.trees] # choose the prediction given by the maximum number of trees used   

    cb_rloss_res <- rloss(tau.pred = tau_cb_pred, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(cb_rloss_res$nnls_coeff, cb_rloss_res$mse, cb_rloss_res$mse_debiased))


    #
    # CausalMARS
    #

    cm <- causalMARS(x = X, tx = W, y = Y)
    tau_cm_pred <- predict(cm, X)
    cm_rloss_pred <- rloss(tau.pred = tau_cm_pred, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(cm_rloss_pred$nnls_coeff, cm_rloss_pred$mse, cm_rloss_pred$mse_debiased))


    #
    # PTOForest
    #
    ptof <- PTOforest(x = X, tx = W, y = Y)
    tau_ptof_pred <- predict(ptof, X)
    ptof_rloss_pred <- rloss(tau.pred = tau_ptof_pred, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(ptof_rloss_pred$nnls_coeff, ptof_rloss_pred$mse, ptof_rloss_pred$mse_debiased))



    ##################################################
    # X-learner
    ##################################################
    xb <- xboost(x = X, w = W, y = Y)
    tau_xb_pred <- predict(xb, X)

    xb_rloss <- rloss(tau.pred = tau_xb_pred, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(xb_rloss$nnls_coeff, xb_rloss$mse, xb_rloss$mse_debiased))


    x_rf <- X_RF(feat = X, tr = W, yobs = Y)
    x_rf_tau <- EstimateCate(x_rf, X)
    x_rf_rloss <- rloss(tau.pred = x_rf_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(x_rf_rloss$nnls_coeff, x_rf_rloss$mse, x_rf_rloss$mse_debiased))


    x_bart <- X_BART(feat = X, tr = W, yobs = Y)
    x_bart_tau <- EstimateCate(x_bart, X)
    x_bart_rloss <- rloss(tau.pred = x_bart_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(x_bart_rloss$nnls_coeff, x_bart_rloss$mse, x_bart_rloss$mse_debiased))


    ##################################################
    # R-learner
    ##################################################
    rb <- rboost(X, W, Y, p_hat = prob, m_hat = Y.hat, nthread = 4)
    rb_tau <- predict(rb)

    rl <- rlasso(X, W, Y, p_hat = prob, m_hat = Y.hat)
    rl_tau <- predict(rl)

    print(round(c(mean((Y - Y.hat - rb_tau * (W - prob))^2), mean((Y - Y.hat - rl_tau * (W - prob))^2)), 4)) # boosting wins

    rb_rloss <- rloss(tau.pred = rb_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(rb_rloss$nnls_coeff, rb_rloss$mse, rb_rloss$mse_debiased))

    # RS learner
    rl_RS <- rlasso(X, W, Y, p_hat = prob, m_hat = Y.hat, rs = TRUE) # returned same predictions for all patients
    rl_RS_tau <- predict(rl_RS)[,1]
    rs_rloss <- rloss(tau.pred = rl_RS_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    # do not include

    ##################################################
    # R-stack
    ##################################################

    RESP = Y - Y.hat
    R.mat <- cbind(1, W - prob,
                    (W - prob) * cf_tau_pred$predictions,
                    (W - prob) * tau_cb_pred,
                    (W - prob) * tau_cm_pred,
                    (W - prob) * tau_ptof_pred,
                    (W - prob) * tau_xb_pred,
                    (W - prob) * x_bart_tau,
                    (W - prob) * x_rf_tau,
                    (W - prob) * rb_tau
                    )

    stack = nnls(as.matrix(R.mat), RESP, constrained = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
    print("coefs")
    print(stack)

    tau.stack = stack[2] +
        stack[3] * cf_tau_pred$predictions +
        stack[4] * tau_cb_pred +
        stack[5] * tau_cm_pred +
        stack[6] * tau_ptof_pred +
        stack[7] * tau_xb_pred + 
        stack[8] * x_bart_tau +
        stack[9] * x_rf_tau +
        stack[10] * rb_tau

    mean((Y - Y.hat - tau.stack * (W - prob))^2)

    tau_stack_rloss <- rloss(tau.pred = tau.stack, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(tau_stack_rloss$nnls_coeff, tau_stack_rloss$mse, tau_stack_rloss$mse_debiased))


    # reformat comparison table
    method_list <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XBART", "XRF", "RBoosting", "stack")
    rownames(compare_rloss) <- method_list
    colnames(compare_rloss) <- c("b", "c", "alpha", "mse", "mse_debiased")

    write.csv(compare_rloss, "./res/NCT00364013_model_sel_res.csv", row.names = TRUE)

}