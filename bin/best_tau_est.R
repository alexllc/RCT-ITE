source("./bin/load_lib.R")
source("./bin/hyperparameter_tuning.R")
source("./bin/mod_PTOforest.R")
source("./bin/mod_causalboost.R")

prob = 0.5
Q = 4
tuned_cb_param = TRUE
tuned_cm_param = TRUE
tune_ptof_param = TRUE
perform_xb = FALSE

# Function to find the best tau predictor using R-loss criteria averaging over several folds
find_best_tau_estimator <- funciton(Y = NULL, X = NULL, W = NULL, prob = 0.5, Q = 4, tuned_cb_param = TRUE, perform_xb = FALSE) {

    # For testing only:
    # Done 
    # NCT00339183, NCT00364013

    trial <- "NCT00113763"
    source(paste0("./bin/load_CRC_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    Y <- get(trial)[[2]][[1]]
    W <- get(trial)[[3]]


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
    while (var(cf_tau_pred$predictions) == 0) { # tuned forests sometimes give equal predictions, if that's the case, rebuild forest
        cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
        cf_tau_pred <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)
    }

    cf_mse <- mean((ho_Y - Y.hat - cf_tau_pred$predictions * (ho_W - prob))^2)
    cf_rloss_res <- rloss(tau.pred = cf_tau_pred$predictions, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)

    compare_rloss <- rbind(compare_rloss, c(cf_rloss_res$nnls_coeff, cf_rloss_res$mse, cf_rloss_res$mse_debiased))


    ##################################################
    # Powers: CausalBoost, CausalMARS, PTOForest
    ##################################################

    if(tuned_cb_param) {
        cb_param <- read.csv("./dat/NCT00339183_cb_param.csv")
        # cb_param <- find_cb_param(train_X, train_Y, train_W, num_search_rounds = 10)
        cb_param <- cb_param[which(cb_param$mean_cvm_effect == min(cb_param$mean_cvm_effect)),]
        cb_param <- list(num.trees = cb_param$num_trees_min_effect, maxleaves = cb_param$max_leaves, eps = cb_param$eps, splitSpread = cb_param$split_spread)
    } else { # use default setting
        cb_param <- list(num.trees = 500, maxleaves = 4, eps = 0.01, splitSpread = 0.1)
    }
    # now we estimate tau from CausalBoosting
    
    cb <- do.call(causalBoosting, append(list(x = train_X, tx = train_W, y = train_Y), as.list(cb_param)))
    cb_tau_pred <- predict(cb, newx = ho_X)
    cb_tau_pred <- cb_tau_pred[,cb_param$num.trees] # choose the prediction given by the maximum number of trees used   

    cb_rloss_res <- rloss(tau.pred = cb_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(cb_rloss_res$nnls_coeff, cb_rloss_res$mse, cb_rloss_res$mse_debiased))


    #
    # CausalMARS
    #
    if (tuned_cm_param){
        cm_cv_performance <- find_cm_param(x = train_X, tx = train_W, y = train_Y, verbose = TRUE)
        cm_param_used <- as.list(cm_cv_performance[1,1:3]) # choose smallest first
        cm <- try(do.call(causalMARS, append(list(x = train_X, tx = train_W, y = train_Y), cm_param_used)))
        tau_cm_pred <- try(predict(cm, ho_X))

        param_counter <- 1
        while( class(cm) == "try-error" | class(tau_cm_pred) == "try-error") { # sometimes you still get singular matricies for the tuned parameters, so tune it until it can fit a CM model. The param_counter determines when you should stop trying
            cm_param_used <- as.list(cm_param[param_counter + 1, 1:3]) # choose the next best param set
            cm <- try(do.call(causalMARS, append(list(x = train_X, tx = train_W, y = train_Y), cm_param_used)))
            tau_cm_pred <- try(predict(cm, ho_X))
            param_counter <- param_counter + 1

            # re-evaluate parameters if too many of the tested parameters failed
            if(param_counter > ceiling(dim(cm_param)[1] * 0.1)) { 
                cm_param <- find_cm_param(x = train_X, tx = train_W, y = train_Y, verbose = TRUE)
                param_counter <- 0
            }
        }

    } else {
        cm <- causalMARS(x = X, tx = W, y = Y)
    }

    cm_rloss_pred <- rloss(tau.pred = tau_cm_pred, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(cm_rloss_pred$nnls_coeff, cm_rloss_pred$mse, cm_rloss_pred$mse_debiased))


    #
    # PTOForest
    #
    if (tune_ptof_param) {
        ptof_param_search <- find_ptof_param(x = train_X, y = train_Y, tx = train_W, validation_fold = 4, num_search_rounds = 50)
        ptof_param <- ptof_param_search[1, 1:3]
        ptof_param <- list(num.trees = ptof_param$num_trees, mtry = ptof_param$mtry, min.node.size = ptof_param$min_node_size)
        ptof <- do.call(PTOforest, append(list(x = train_X, y = train_Y, tx = train_W, postprocess = TRUE, verbose = TRUE), ptof_param))
    } else {
        ptof <- PTOforest(x = train_X, tx = train_W, y = train_Y, verbose = TRUE)

    }
    ptof_tau_pred <- predict(ptof, ho_X)
    ptof_rloss_pred <- rloss(tau.pred = ptof_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(ptof_rloss_pred$nnls_coeff, ptof_rloss_pred$mse, ptof_rloss_pred$mse_debiased))



    ##################################################
    # X-learner
    ##################################################
    if (perform_xb) {
        xb <- xboost(x = train_X, w = train_W, y = train_Y, verbose = TRUE)
        xb_tau_pred <- predict(xb, ho_X)
        xb_rloss <- rloss(tau.pred = xb_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
        compare_rloss <- rbind(compare_rloss, c(xb_rloss$nnls_coeff, xb_rloss$mse, xb_rloss$mse_debiased))
    } else {
        xb_tau_pred <- rep(0, dim(ho_X)[1])
        compare_rloss <- rbind(compare_rloss, rep(NA, 3))
    }


    x_rf <- X_RF(feat = train_X, tr = train_W, yobs = train_Y)
    x_rf_tau <- EstimateCate(x_rf, ho_X)
    x_rf_rloss <- rloss(tau.pred = x_rf_tau, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(x_rf_rloss$nnls_coeff, x_rf_rloss$mse, x_rf_rloss$mse_debiased))


    x_bart <- X_BART(feat = train_X, tr = train_W, yobs = train_Y)
    x_bart_tau <- EstimateCate(x_bart, ho_X)
    x_bart_rloss <- rloss(tau.pred = x_bart_tau, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(x_bart_rloss$nnls_coeff, x_bart_rloss$mse, x_bart_rloss$mse_debiased))


    ##################################################
    # R-learner
    ##################################################

    train_Y_hat_boost <- predict(Y.boost)

    train_Y_hat_lasso <- predict(Y.lasso, newx = train_X)[,1]

    # which method has a smaller CV error?
    print(round(c(RMSE(train_Y_hat_boost, train_Y), RMSE(train_Y_hat_lasso, train_Y)), 4))
    
    if( RMSE(train_Y_hat_boost, train_Y) < RMSE(train_Y_hat_boost, train_Y) ) {
        train_Y_hat <- train_Y_hat_boost
    } else {
        train_Y_hat <- train_Y_hat_lasso
    }

    rb <- rboost(train_X, train_W, train_Y, p_hat = prob, m_hat = train_Y_hat, nthread = 40)
    rb_tau <- predict(rb, ho_X)
    rb_rloss <- rloss(tau.pred = rb_tau, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(rb_rloss$nnls_coeff, rb_rloss$mse, rb_rloss$mse_debiased))

    rl <- rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = train_Y_hat)
    rl_tau <- predict(rl, ho_X)[,1]
    rl_rloss <- try(rloss(tau.pred = rl_tau, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob))
    if (class(rl_rloss) == "try-error") {
        rl_tau <- rep(0, dim(ho_X)[1])
        compare_rloss <- rbind(compare_rloss, rep(NA, 3))
    } else {
        compare_rloss <- rbind(compare_rloss, c(rl_rloss$nnls_coeff, rl_rloss$mse, rl_rloss$mse_debiased))
    }

    print(round(c(mean((train_Y - train_Y_hat - rb_tau * (train_W - prob))^2), mean((train_Y - train_Y_hat - rl_tau * (train_W - prob))^2)), 4)) # boosting wins


    # RS learner
    rl_RS <- rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = train_Y_hat, rs = TRUE) # returned same predictions for all patients
    rl_RS_tau <- predict(rl_RS, ho_X)[,1]
    rs_rloss <- try(rloss(tau.pred = rl_RS_tau, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob))
    if (class(rs_rloss) == "try-error") {
        rl_RS_tau <- rep(0, dim(ho_X)[1])
        compare_rloss <- rbind(compare_rloss, rep(NA, 3))
    } else {
        compare_rloss <- rbind(compare_rloss, rep(NA, 3))
    }

    ##################################################
    # R-stack
    ##################################################

    RESP = train_Y - Y.hat
    R.mat <- cbind(1, train_W - prob,
                    (train_W - prob) * cf_tau_pred$predictions,
                    (train_W - prob) * cb_tau_pred,
                    (train_W - prob) * tau_cm_pred,
                    (train_W - prob) * ptof_tau_pred,
                    # (train_W - prob) * xb_tau_pred,
                    (train_W - prob) * x_bart_tau,
                    (train_W - prob) * x_rf_tau,
                    (train_W - prob) * rb_tau
                    # (train_W - prob) * rl_tau
                    # (train_W - prob) * rl_RS_tau
                    )

    stack = nnls(as.matrix(R.mat), RESP, constrained = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,  TRUE))
    print("coefs")
    print(stack)

    tau.stack = stack[2] +
        stack[3] * cf_tau_pred$predictions +
        stack[4] * cb_tau_pred +
        stack[5] * tau_cm_pred +
        stack[6] * ptof_tau_pred +
        # stack[7] * tau_xb_pred + 
        stack[7] * x_bart_tau +
        stack[8] * x_rf_tau +
        stack[9] * rb_tau

    mean((train_Y - train_Y_hat - tau.stack * (train_W - prob))^2)

    tau_stack_rloss <- rloss(tau.pred = tau.stack, Y = ho_Y, W = ho_W, Y.hat = Y.hat, prob = prob)
    compare_rloss <- rbind(compare_rloss, c(tau_stack_rloss$nnls_coeff, tau_stack_rloss$mse, tau_stack_rloss$mse_debiased))


    # reformat comparison table
    method_list <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XRF", "XBART", "RBoosting", "RLasso", "RSlearner", "stack")
    rownames(compare_rloss) <- method_list
    colnames(compare_rloss) <- c("b", "c", "alpha", "mse", "mse_debiased")

    write.csv(compare_rloss, paste0("./res/", trial, "_model_sel_res.csv"), row.names = TRUE)

}