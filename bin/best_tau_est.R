# Using R-loss criteria to find the best tau estimator 2022-04-13

source("./bin/load_lib.R")
source("./bin/hyperparameter_tuning.R")
source("./bin/mod_PTOforest.R")
source("./bin/mod_causalboost.R")

prob = 0.5
Q = 4
tuned_cb_param = FALSE
tuned_cm_param = FALSE
tune_ptof_param = FALSE
perform_xb = FALSE
# 
# trial_ls <- c("NCT00364013", "NCT00339183", "NCT00115765", "NCT00113763", "NCT00079274",
#                 "NCT00460265",
#                 "NCT00041119_length", "NCT00041119_chemo",
#                 "NCT00003299", "NCT00119613")

# done: "NCT00052910", "NCT00113763", "NCT00460265","NCT00364013",  "NCT00115765", "NCT00339183"
# pending: "NCT00115765_oxa", 
trial_ls <- c("NCT00041119_length")

# Function to find the best tau predictor using R-loss criteria averaging over several folds
# find_best_tau_estimator <- funciton(Y = NULL, X = NULL, W = NULL, prob = 0.5, Q = 4, tuned_cb_param = TRUE, tuned_cm_param = TRUE, tune_ptof_param = TRUE, perform_xb = TRUE) {

for (trial in trial_ls) {

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
    # outcome_list <- get(paste0(trial, "_outcomes"))
    outcome_list <- c("OS", "RFS")

    for (outcome in outcome_list) {

        message(paste0("Processing outcome: ", outcome))

        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }

        # We need to take out 1/Q sample as the holdout test data, the rest of the samples are used for training
        
        ho_matrix <- read.csv(paste0("./dat/cv_ho_ids/est_tau_", trial, "_holdout_id.csv"))

        for (iter in 1:Q) {
            
            tau_ls <- list()
            compare_rloss <- data.frame(b = numeric(), c = numeric(), alpha = numeric(), mse = numeric(), debiased_mse = numeric())

            message(paste0(rep("=", 80)))
            message(paste0("CV iteration number: ", iter))
            message(paste0(rep("=", 80)))
            
            holdout_sample <- ho_matrix[,iter]
            available_sample <- seq(1:dim(X)[1])
            train_sample <- available_sample[!available_sample %in% holdout_sample]

            ho_X <- X[holdout_sample,]
            ho_W <- W[holdout_sample]
            ho_Y <- Y[holdout_sample]

            train_X <- X[train_sample,]
            train_W <- W[train_sample]
            train_Y <- Y[train_sample]

            # Since CausalBoost is not implemented w.r.t. the causal forest in Imbens and Athey, there is not debiased.error to take out of context

            # let's first fit cvlasso or cv boosting to estimate nuisance component marginal response model (m.hat)

            message(paste0(rep("=", 80)))
            message(paste0("Estimating m hat for holdout data."))
            message(paste0(rep("=", 80)))

            Y_boost <- cvboost(train_X, train_Y, objective = "reg:squarederror", nthread = 40) # non-binary outcome
            Y_hat_boost <- predict(Y_boost, newx = ho_X)

            Y_lasso <- cv.glmnet(train_X, train_Y, keep = TRUE, family = "gaussian")
            Y_hat_lasso <- predict(Y_lasso, newx = ho_X)[,1]

            # which method has a smaller CV error?
            print("RMSE of m hat estimated by boosting vs LASSO:")
            print(round(c(RMSE(Y_hat_boost, ho_Y), RMSE(Y_hat_lasso, ho_Y)), 4))
            
            if( RMSE(Y_hat_boost, ho_Y) < RMSE(Y_hat_lasso, ho_Y) ) {
                Y_hat <- Y_hat_boost
            } else {
                Y_hat <- Y_hat_lasso
            }
            #
            # Causal forest
            #

            # R-loss already implemented by the package and is contained in the debiased.error column in predict.causal_forest(). You can obtain MSE by mean(cf$debiased.error), the values were already squared.
            # since causal trees are honest already, there's no need to further split into k-fold estimation

            message(paste0(rep("=", 80)))
            message(paste0("Building Causal Forest."))
            message(paste0(rep("=", 80)))

            cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
            cf_tau_pred <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)
            while (var(cf_tau_pred$predictions) == 0) { # tuned forests sometimes give equal predictions, if that's the case, rebuild forest
                cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
                cf_tau_pred <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)
            }

            tau_ls <- append(tau_ls, list(cf_tau_pred$predictions))

            cf_mse <- mean((ho_Y - Y_hat - cf_tau_pred$predictions * (ho_W - prob))^2)
            cf_rloss_res <- rloss(tau.pred = cf_tau_pred$predictions, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)

            compare_rloss <- rbind(compare_rloss, c(cf_rloss_res$nnls_coeff, cf_rloss_res$mse, cf_rloss_res$mse_debiased))


            ##################################################
            # Powers: CausalBoost, CausalMARS, PTOForest
            ##################################################

            message(paste0(rep("=", 80)))
            message(paste0("Performing causal boosting."))
            message(paste0(rep("=", 80)))

            if(tuned_cb_param) {
                cb_param <- find_cb_param(train_X, train_Y, train_W, num_search_rounds = 10)
                cb_param <- cb_param[which(cb_param$mean_cvm_effect == min(cb_param$mean_cvm_effect)),]
                cb_param <- list(num.trees = cb_param$num_trees_min_effect, maxleaves = cb_param$max_leaves, eps = cb_param$eps, splitSpread = cb_param$split_spread)
            } else { # use default setting
                cb_param_tune_res <- read.csv(paste0("./dat/cb_param/", trial, "_", imp_type_Y, "_", outcome, "_cb_param.csv")) # if skipping parameters tuning, we will use the previously tuned parameter for this trial
                cb_param_tune_res <- filter(cb_param_tune_res, num_trees_min_effect != 1) # avoid using parameteres with only 1 boosting tree
                cb_param <- cb_param_tune_res[which.min(cb_param_tune_res$mean_cvm_effect),]
                cb_param <- list(num.trees = cb_param$num_trees_min_effect, maxleaves = cb_param$max_leaves, eps = cb_param$eps, splitSpread = cb_param$split_spread)
            }
            # now we estimate tau from CausalBoosting
            
            cb <- do.call(causalBoosting, append(list(x = train_X, tx = train_W, y = train_Y), as.list(cb_param)))
            cb_tau_pred <- predict(cb, newx = ho_X)
            cb_tau_pred <- cb_tau_pred[,cb_param$num.trees] # choose the prediction given by the maximum number of trees used   

            tau_ls <- append(tau_ls, list(cb_tau_pred))

            cb_rloss_res <- rloss(tau.pred = cb_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
            compare_rloss <- rbind(compare_rloss, c(cb_rloss_res$nnls_coeff, cb_rloss_res$mse, cb_rloss_res$mse_debiased))


            #
            # CausalMARS
            #
            message(paste0(rep("=", 80)))
            message(paste0("Building Performing causal MARS."))
            message(paste0(rep("=", 80)))
            if (tuned_cm_param | !file.exists(paste0("./res/params/", trial, "_cm_params.csv"))){
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
                write.csv(cm_cv_performance, paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_cm_params.csv"), row.names = FALSE)
            } else {
                cm_cv_performance <- read.csv(paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_cm_params.csv"))
                cm_param_used <- as.list(cm_cv_performance[1,1:3]) # choose smallest first
                cm <- do.call(causalMARS, append(list(x = train_X, tx = train_W, y = train_Y), cm_param_used))
                tau_cm_pred <- try(predict(cm, ho_X))
            }

            tau_ls <- append(tau_ls, list(tau_cm_pred))

            cm_rloss_pred <- rloss(tau.pred = tau_cm_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
            compare_rloss <- rbind(compare_rloss, c(cm_rloss_pred$nnls_coeff, cm_rloss_pred$mse, cm_rloss_pred$mse_debiased))


            #
            # PTOForest
            #
            message(paste0(rep("=", 80)))
            message(paste0("Building PTOForest."))
            message(paste0(rep("=", 80)))
            if (tune_ptof_param | !file.exists(paste0("./res/params/", trial, "_ptof_params.csv"))) {
                ptof_param_search <- find_ptof_param(x = train_X, y = train_Y, tx = train_W, validation_fold = 4, num_search_rounds = 50)
                ptof_param <- ptof_param_search[1, 1:3]
                ptof_param <- list(num.trees = ptof_param$num_trees, mtry = ptof_param$mtry, min.node.size = ptof_param$min_node_size)
                ptof <- try(do.call(PTOforest, append(list(x = train_X, y = train_Y, tx = train_W, postprocess = FALSE, verbose = TRUE), ptof_param))) # sometimes PTO forest can fail to produce some estimates with holdout data
                write.csv(ptof_param_search, paste0("./res/params/", trial, "_ptof_params.csv"), row.names = FALSE)
            } else {
                ptof_param_search <- read.csv(paste0("./res/params/", trial, "_ptof_params.csv"))
                ptof_param <- ptof_param_search[1, 1:3]
                ptof_param <- list(num.trees = ptof_param$num_trees, mtry = ptof_param$mtry, min.node.size = ptof_param$min_node_size)
                ptof <- try(do.call(PTOforest, append(list(x = train_X, y = train_Y, tx = train_W, postprocess = FALSE, verbose = TRUE), ptof_param))) # sometimes PTO forest can fail to produce some estimates with holdout data
            }
            ptof_tau_pred <- try(predict(ptof, ho_X))
            ptof_rloss_pred <- try(rloss(tau.pred = ptof_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob))

            if (class(ptof) == "try-error" | class(ptof_tau_pred) == "try-error" | class(ptof_rloss_pred) == "try-error") {
                print("Cannot build this model.")
                ptof_tau_pred <- rep(0, dim(ho_X)[1])
                compare_rloss <- rbind(compare_rloss, rep(NA, 3))
            } else {
                compare_rloss <- rbind(compare_rloss, c(ptof_rloss_pred$nnls_coeff, ptof_rloss_pred$mse, ptof_rloss_pred$mse_debiased))
            }

            tau_ls <- append(tau_ls, list(ptof_tau_pred))




            ##################################################
            # X-learner
            ##################################################
            #
            # X-learner of boosting
            #
            message(paste0(rep("=", 80)))
            message(paste0("Building X-learner."))
            message(paste0(rep("=", 80)))

            if (perform_xb) {
                xb <- xboost(x = train_X, w = train_W, y = train_Y, ntrees_max= 200, num_search_rounds = 5, verbose = TRUE)
                xb_tau_pred <- predict(xb, ho_X)
                xb_rloss <- rloss(tau.pred = xb_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
                compare_rloss <- rbind(compare_rloss, c(xb_rloss$nnls_coeff, xb_rloss$mse, xb_rloss$mse_debiased))
            } else {
                xb_tau_pred <- rep(0, dim(ho_X)[1])
                compare_rloss <- rbind(compare_rloss, rep(NA, 3))
            }

            tau_ls <- append(tau_ls, list(xb_tau_pred))
            
            #
            # X-learner of random forests
            #

            x_rf <- X_RF(feat = train_X, tr = train_W, yobs = train_Y)
            x_rf_tau <- EstimateCate(x_rf, ho_X)
            x_rf_rloss <- rloss(tau.pred = x_rf_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
            compare_rloss <- rbind(compare_rloss, c(x_rf_rloss$nnls_coeff, x_rf_rloss$mse, x_rf_rloss$mse_debiased))

            tau_ls <- append(tau_ls, list(x_rf_tau))

            #
            # X-learner of BART
            #

            x_bart <- X_BART(feat = train_X, tr = train_W, yobs = train_Y)
            x_bart_tau <- EstimateCate(x_bart, ho_X)
            x_bart_rloss <- rloss(tau.pred = x_bart_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
            compare_rloss <- rbind(compare_rloss, c(x_bart_rloss$nnls_coeff, x_bart_rloss$mse, x_bart_rloss$mse_debiased))

            tau_ls <- append(tau_ls, list(x_bart_tau))

            ##################################################
            # R-learner
            ##################################################
            message(paste0(rep("=", 80)))
            message(paste0("Building R-learner."))
            message(paste0(rep("=", 80)))

            train_Y_hat_boost <- predict(Y_boost)

            train_Y_hat_lasso <- predict(Y_lasso, newx = train_X)[,1]

            # which method has a smaller CV error?
            print(round(c(RMSE(train_Y_hat_boost, train_Y), RMSE(train_Y_hat_lasso, train_Y)), 4))
            
            if( RMSE(train_Y_hat_boost, train_Y) < RMSE(train_Y_hat_boost, train_Y) ) {
                train_Y_hat <- train_Y_hat_boost
            } else {
                train_Y_hat <- train_Y_hat_lasso
            }

            #
            # R-learner of boosting
            #

            rb <- rboost(train_X, train_W, train_Y, p_hat = prob, m_hat = train_Y_hat, nthread = 40)
            rb_tau <- predict(rb, ho_X)
            rb_rloss <- rloss(tau.pred = rb_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
            compare_rloss <- rbind(compare_rloss, c(rb_rloss$nnls_coeff, rb_rloss$mse, rb_rloss$mse_debiased))

            tau_ls <- append(tau_ls, list(rb_tau))

            #
            # R-learner of LASSO
            #

            rl <- rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = train_Y_hat)
            rl_tau <- predict(rl, ho_X)[,1]
            rl_rloss <- try(rloss(tau.pred = rl_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob))
            if (class(rl_rloss) == "try-error") {
                print("Cannot build this model.")
                rl_tau <- rep(0, dim(ho_X)[1])
                compare_rloss <- rbind(compare_rloss, rep(NA, 3))
            } else {
                compare_rloss <- rbind(compare_rloss, c(rl_rloss$nnls_coeff, rl_rloss$mse, rl_rloss$mse_debiased))
            }
            tau_ls <- append(tau_ls, list(rl_tau))
            print(round(c(mean((train_Y - train_Y_hat - rb_tau * (train_W - prob))^2), mean((train_Y - train_Y_hat - rl_tau * (train_W - prob))^2)), 4)) # boosting wins

            #
            # RS-learner of LASSO
            #

            rl_RS <- rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = train_Y_hat, rs = TRUE) # returned same predictions for all patients
            rl_RS_tau <- predict(rl_RS, ho_X)[,1]
            rs_rloss <- try(rloss(tau.pred = rl_RS_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob))
            if (class(rs_rloss) == "try-error") {
                print("Cannot build this model.")
                rl_RS_tau <- rep(0, dim(ho_X)[1])
                compare_rloss <- rbind(compare_rloss, rep(NA, 3))
            } else {
                compare_rloss <- rbind(compare_rloss, c(rs_rloss$nnls_coeff, rs_rloss$mse, rs_rloss$mse_debiased))
            }
            tau_ls <- append(tau_ls, list(rl_RS_tau))

            names(tau_ls) <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XRF", "XBART", "RBoosting", "RLasso", "RSlearner")
            
            ##################################################
            # R-stack
            ##################################################

            RESP <- ho_Y - Y_hat
            R_mat <- cbind(1, ho_W - rep(prob, length(ho_W)))
            tau_counter <- 0
            for(tau_pred in tau_ls) {
                if (sum(tau_pred) != 0){
                    R_mat <- cbind(R_mat, (ho_W - rep(prob, length(ho_W))) * tau_pred)
                    tau_counter <- tau_counter + 1
                }
            }

            stack <- try(nnls(as.matrix(R_mat), RESP, constrained = c(FALSE, FALSE, rep(TRUE, tau_counter))))

            if (class(stack) == "try-error") {
                print("Cannot build this model.")
                tau_stack <- rep(0, dim(ho_X)[1])
                compare_rloss <- rbind(compare_rloss, rep(NA, 3))
            } else {
                print("coefs")
                print(stack)

                tau_stack <- stack[2]
                j <- 3
                for (i in 1:length(tau_ls)) {
                    if (sum(tau_ls[[i]]) != 0){
                        tau_stack <- tau_stack + stack[j] * tau_ls[[i]]
                        j <- j + 1
                    }
                }

                print("R-loss of the R-stack learner: ")
                mean((ho_Y - Y_hat - tau_stack * (ho_W - rep(prob, length(ho_W))))^2)

                tau_stack_rloss <- try(rloss(tau.pred = tau_stack, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob))
                if(class(tau_stack_rloss) == "try-error") {
                    rstack_tau <- rep(0, dim(ho_X)[1])
                    compare_rloss <- rbind(compare_rloss, rep(NA, 3))
                } else {
                    compare_rloss <- rbind(compare_rloss, c(tau_stack_rloss$nnls_coeff, tau_stack_rloss$mse, tau_stack_rloss$mse_debiased))
                }
            }


            # reformat comparison table
            method_list <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XRF", "XBART", "RBoosting", "RLasso", "RSlearner", "stack")
            rownames(compare_rloss) <- method_list
            colnames(compare_rloss) <- c("b", "c", "alpha", "mse", "mse_debiased")

            write.csv(compare_rloss, paste0("./res/", trial, "_model_sel_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"), row.names = TRUE)
        } # end of CV fold loop
    } # end of outcome type loop
} # end of trial loop