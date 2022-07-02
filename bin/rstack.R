#' Function to calculate tau estimate using out-of-fold cross-fitting described in R-stack in Nie & Wager (2021).

source("./bin/load_lib.R")
source("./bin/hyperparameter_tuning.R")
source("./bin/mod_PTOforest.R")
source("./bin/mod_causalboost.R")

#' The methods considered are : causal forest, causal boosting, causal MARS, PTO Forest, X-learners (boosting, random forest, BART), R-learners (boosting, LASSO, RS learner)
#' The tau estimates are all fitted using a model built from out-of-fold training sample set. Each of the fold takes turn to be the testing set, so this R-stack is performed Q times depending on the number of cross-fitting folds.
#' All methods concerned are tuned according to CV R-loss values. Causal boosting takes a long time to tune, so it is preferable to load a previously tuned set of parameteres (and corresponding training vs testing sample set) rather than tuning it everytime.
#' Some of the tau estimators may be invalid for the dataset (e.g. the covariate matrix is not singular), in this case the R-stack will exclude the said method entirely.
#' @param Y numeric vector of outcome
#' @param X numeric matrix of covariates without missing values
#' @param W numeric vector of binary treatment, 1: treated, 0: control
#' @param prob numeric vector or value indicating propensity score
#' @param Q numeric value number of folds
#' @param tuned_cb_param bool whether casual boosting parameters should be tuned or loaded
#' @param tuned_cm_param bool whether causal MARS parameters should be tuned or use default
#' @param tune_ptof_param bool whether PTO forest parameters should be tuned or use default 
#' @param perform_xb bool whether X-learner boosting should be performed. This is skipped in our analysis because it requires huge resources and time without the guarantee of good performance.

perform_rstack <- function(Y = NULL, X = NULL, W = NULL, trial_name = NULL, outcome = NULL, imp_type_Y = NULL, prob = 0.5, Q = 4, tuned_cb_param = FALSE, tuned_cm_param = TRUE, tune_ptof_param = TRUE, perform_xb = FALSE, perform_cb = FALSE) {

    available_sample <- 1:dim(X)[1]
    cf_tau <- data.frame(ID = available_sample)
    rownames(cf_tau) <- cf_tau$ID

    # Load previously generated training and testing IDs for reprodcibility
    ho_matrix <- read.csv(paste0("./dat/cv_ho_ids/est_tau_", trial, "_holdout_id.csv"))

    for (iter in 1:Q) {

        message(paste0(rep("=", 80)))
        message(paste0("CV iteration number: ", iter))
        message(paste0(rep("=", 80)))
        
        tau_ls <- list()
        
        compare_rloss <- data.frame(b = numeric(), c = numeric(), alpha = numeric(), mse = numeric(), debiased_mse = numeric())

        ho_id <- ho_matrix[,iter]
        available_sample <- seq(1:dim(X)[1])
        # holdout_sample <- sample(seq(1:dim(X)[1]), floor(dim(X)[1] / 4))
        train_id <- available_sample[!available_sample %in% ho_id] 

        train_X <- X[train_id,]
        train_Y <- Y[train_id]
        train_W <- W[train_id]

        ho_X <- X[ho_id,]
        ho_Y <- Y[ho_id]
        ho_W <- W[ho_id]

        #
        # R-stack
        #
        if (!(imp_type_Y == "efronYn" | imp_type_Y == "efron")) {
            binary_Y <- TRUE
        } else {
            binary_Y <- FALSE
        }
        Y_hat <- get_mhat(X_train = train_X, W_train = train_W, Y_train = train_Y, X_ho = ho_X, Y_ho = ho_Y, binary_Y = binary_Y)
        
        #
        # Causal forest
        #
        message(paste0(rep("=", 80)))
        message(paste0("Building Causal Forest."))
        message(paste0(rep("=", 80)))

        cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
        cf_tau <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)
        cf_tau_pred <- cf_tau$predictions
        cf_rep_count <- 0
        while (var(cf_tau_pred) == 0 & cf_rep_count < 5) { # tuned forests sometimes give equal predictions, if that's the case, rebuild forest
            cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
            cf_tau <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)
            cf_tau_pred <- cf_tau$predictions
            cf_rep_count <- cf_rep_count + 1
        }

        if (var(cf_tau_pred$predictions) == 0 & cf_rep_count >= 5) {
            message("Failed to build a tuned CF with varying predictions.")
            cf_tau_pred <- rep(0, dim(ho_X)[1])
            compare_rloss <- rbind(compare_rloss, rep(NA, 3))
        } else {
            cf_mse <- mean((ho_Y - Y_hat - cf_tau_pred * (ho_W - prob))^2)
            cf_rloss_res <- rloss(tau.pred = cf_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)

            compare_rloss <- rbind(compare_rloss, c(cf_rloss_res$nnls_coeff, cf_rloss_res$mse, cf_rloss_res$mse_debiased))

            write.csv(cf_tau_pred, file = paste0("./res/indiv_model_tau/", trial, "_CF_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)
        }

        ##################################################
        # Powers: CausalBoost, CausalMARS, PTOForest
        ##################################################

        if (perform_cb) {
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

            cb_rloss_res <- rloss(tau.pred = cb_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
            compare_rloss <- rbind(compare_rloss, c(cb_rloss_res$nnls_coeff, cb_rloss_res$mse, cb_rloss_res$mse_debiased))
        } else {
            print("Skipping Causal Boosting.")
            cb_tau_pred <- rep(0, dim(ho_X)[1])
            compare_rloss <- rbind(compare_rloss, rep(NA, 3))
        }


        #
        # CausalMARS
        #
        message(paste0(rep("=", 80)))
        message(paste0("Building Performing causal MARS."))
        message(paste0(rep("=", 80)))
        if (tuned_cm_param | !file.exists( paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_", iter,"_cm_params.csv"))){
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
            write.csv(cm_cv_performance, paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_", iter,"_cm_params.csv"), row.names = FALSE)

        } else {
            cm_cv_performance <- read.csv(paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_", iter, "_cm_params.csv"))
            cm_param_used <- as.list(cm_cv_performance[1,1:3]) # choose smallest first
            cm <- do.call(causalMARS, append(list(x = train_X, tx = train_W, y = train_Y), cm_param_used))
            tau_cm_pred <- try(predict(cm, ho_X))
        }

        cm_rloss_pred <- rloss(tau.pred = tau_cm_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
        compare_rloss <- rbind(compare_rloss, c(cm_rloss_pred$nnls_coeff, cm_rloss_pred$mse, cm_rloss_pred$mse_debiased))

        write.csv(tau_cm_pred, file = paste0("./res/indiv_model_tau/", trial, "_CMARS_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)

        #
        # PTOForest
        #
        message(paste0(rep("=", 80)))
        message(paste0("Building PTOForest."))
        message(paste0(rep("=", 80)))
        rebuild_ptof <- TRUE
        while (rebuild_ptof) {
            if (tune_ptof_param) {
                ptof_param_search <- find_ptof_param(x = train_X, y = train_Y, tx = train_W, validation_fold = 4, num_search_rounds = 50)
                ptof_param <- ptof_param_search[1, 1:3]
                ptof_param <- list(num.trees = ptof_param$num_trees, mtry = ptof_param$mtry, min.node.size = ptof_param$min_node_size)
                ptof <- try(do.call(PTOforest, append(list(x = train_X, y = train_Y, tx = train_W, postprocess = FALSE, verbose = TRUE), ptof_param))) # sometimes PTO forest can fail to produce some estimates with holdout data
            } else {
                ptof <- try(PTOforest(x = train_X, tx = train_W, y = train_Y, verbose = TRUE))
            }

            ptof_tau_pred <- try(predict(ptof, ho_X))
            ptof_rloss_pred <- try(rloss(tau.pred = ptof_tau_pred, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob))

            if (any(is.na(ptof_tau_pred))) {
                message("Out of bag NA predictions generated, will be refitting the model.")
                rebuild_ptof <- any(is.na(ptof_tau_pred))
            } else {
                rebuild_ptof <- FALSE
            }
        }
        write.csv(ptof_param_search, file = paste0(paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_", iter, "_ptof_params_used.csv")))

        if (class(ptof) == "try-error" | class(ptof_tau_pred) == "try-error") {
            print("Cannot build this model.")
            ptof_tau_pred <- rep(0, dim(ho_X)[1])
            compare_rloss <- rbind(compare_rloss, rep(NA, 3))
        } else {
            compare_rloss <- rbind(compare_rloss, c(ptof_rloss_pred$nnls_coeff, ptof_rloss_pred$mse, ptof_rloss_pred$mse_debiased))
        }

        write.csv(ptof_tau_pred, file = paste0("./res/indiv_model_tau/", trial, "_PTOForest_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)


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
        
        
        #
        # X-learner of random forests
        #

        x_rf <- X_RF(feat = train_X, tr = train_W, yobs = train_Y)
        x_rf_tau <- EstimateCate(x_rf, ho_X)

        x_rf_rloss <- rloss(tau.pred = x_rf_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
        compare_rloss <- rbind(compare_rloss, c(x_rf_rloss$nnls_coeff, x_rf_rloss$mse, x_rf_rloss$mse_debiased))

        write.csv(x_rf_tau, file = paste0("./res/indiv_model_tau/", trial, "_XRF_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)

        #
        # X-learner of BART
        #

        x_bart <- X_BART(feat = train_X, tr = train_W, yobs = train_Y)
        x_bart_tau <- EstimateCate(x_bart, ho_X)

        x_bart_rloss <- rloss(tau.pred = x_bart_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
        compare_rloss <- rbind(compare_rloss, c(x_bart_rloss$nnls_coeff, x_bart_rloss$mse, x_bart_rloss$mse_debiased))

        write.csv(x_bart_tau, file = paste0("./res/indiv_model_tau/", trial, "_XBART_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)        


        ##################################################
        # R-learner
        ##################################################
        message(paste0(rep("=", 80)))
        message(paste0("Building R-learner."))
        message(paste0(rep("=", 80)))

        #
        # R-learner of boosting
        #
        
        ## NEED EXTRA M_HAT FOR THESE ALGORITHMS?

        subtrain_Y_hat <- get_mhat(X_train = train_X, W_train = train_W, Y_train = train_Y, X_ho = ho_X, Y_ho = ho_Y, binary_Y = binary_Y, newdat = FALSE)

        rb <- try(rboost(train_X, train_W, train_Y, p_hat = prob, m_hat = subtrain_Y_hat, nthread = 40))
        while(class(rb) == "try-error") {
            rb <- try(rboost(train_X, train_W, train_Y, p_hat = prob, m_hat = Y_hat, nthread = 40))
        }
        rb_tau <- predict(rb, ho_X)

        rb_rloss <- rloss(tau.pred = rb_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob)
        compare_rloss <- rbind(compare_rloss, c(rb_rloss$nnls_coeff, rb_rloss$mse, rb_rloss$mse_debiased))

        write.csv(rb_tau, file = paste0("./res/indiv_model_tau/", trial, "_RBoosting_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)  

        #
        # R-learner of LASSO
        #

        rl <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = subtrain_Y_hat))
        while(class(rl) == "try-error") {
            rl <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = subtrain_Y_hat))
        }
        rl_tau <- predict(rl, ho_X)[,1]

        rl_rloss <- try(rloss(tau.pred = rl_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob))
        if (class(rl_rloss) == "try-error") {
            print("Cannot build this model.")
            rl_tau <- rep(0, dim(ho_X)[1])
            compare_rloss <- rbind(compare_rloss, rep(NA, 3))
        } else {
            compare_rloss <- rbind(compare_rloss, c(rl_rloss$nnls_coeff, rl_rloss$mse, rl_rloss$mse_debiased))
        }
        
        print(round(c(mean((train_Y - subtrain_Y_hat - rb_tau * (train_W - prob))^2), mean((train_Y - subtrain_Y_hat - rl_tau * (train_W - prob))^2)), 4)) 

        write.csv(rl_tau, file = paste0("./res/indiv_model_tau/", trial, "_RLasso_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)  

        #
        # RS-learner of LASSO
        #

        rl_RS <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = subtrain_Y_hat, rs = TRUE)) # returned same predictions for all patients
        while(class(rl_RS) == "try-error") {
            rl_RS <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = subtrain_Y_hat, rs = TRUE))
        }
        rl_RS_tau <- predict(rl_RS, ho_X)[,1]
        rs_rloss <- try(rloss(tau.pred = rl_RS_tau, Y = ho_Y, W = ho_W, Y.hat = Y_hat, prob = prob))
        if (class(rs_rloss) == "try-error") {
            print("Cannot build this model.")
            rl_RS_tau <- rep(0, dim(ho_X)[1])
            compare_rloss <- rbind(compare_rloss, rep(NA, 3))
        } else {
            compare_rloss <- rbind(compare_rloss, c(rs_rloss$nnls_coeff, rs_rloss$mse, rs_rloss$mse_debiased))
        }

        write.csv(rl_RS_tau, file = paste0("./res/indiv_model_tau/", trial, "_RSlearner_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv"), row.names = FALSE)  

        tau_ls <- append(tau_ls, list(cf_tau_pred, cb_tau_pred, tau_cm_pred, ptof_tau_pred, xb_tau_pred, x_rf_tau, x_bart_tau, rb_tau, rl_tau, rl_RS_tau))

        names(tau_ls) <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XRF", "XBART", "RBoosting", "RLasso", "RSlearner")


        ##################################################
        # R-stack
        ##################################################

        RESP <- ho_Y - Y_hat
        R_mat <- cbind(1, ho_W - rep(prob, length(ho_W)))
        tau_counter <- 0
        name_counter <- 1
        tau_matrix <- data.frame(nrow = length(ho_Y))

        for(tau_pred in tau_ls) {
            # nnls will fail due to a non-positive definite matrix if you have a constant tau vector
            # sometimes tree based tau estimators will render NAs in out of bag tau predictions, which will be excluded in the stack
            if (sum(tau_pred, na.rm = TRUE) != 0 & var(tau_pred) != 0 & !any(is.na(tau_pred))){ 
            
                R_mat <- cbind(R_mat, (ho_W - rep(prob, length(ho_W))) * tau_pred)
                tau_counter <- tau_counter + 1

                # Check calibration of ITE estimates
                message(paste0("Calibration of: ", names(tau_ls)[name_counter]))
                calib <- tau_calibration(tau.hat = tau_pred, Y = ho_Y, m.hat = Y_hat, W = ho_W, W.hat = 0.5)
                write.csv(calib, file = paste0("./res/calib/", trial, "_", names(tau_ls)[name_counter], "_calib_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"), row.names = FALSE)
                name_counter <- name_counter + 1

                # Correlation between each methods
                tau_matrix <- cbind(tau_matrix, tau_pred)

            }
        }

        # Save correlations between methods
        tau_matrix$nrow <- NULL
        colnames(tau_matrix) <- names(tau_ls)[sapply(tau_ls, function(X) sum(X) != 0 & var(X) != 0 & !any(is.na(X)))]

        corr_mat <- rcorr(as.matrix(tau_matrix))
        write.csv(corr_mat$P, file = paste0("./res/corr_matrix/", trial, "_corr_matrix_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"), row.names = FALSE)


        stack <- try(nnls(as.matrix(R_mat), RESP, constrained = c(FALSE, FALSE, rep(TRUE, tau_counter))))
        last_stack <- stack

        if (class(stack) == "try-error") {
            print("Cannot build this model, using coeff from the last fold.")
            
            stack <- last_stack
            print(stack)

            tau_stack <- stack[2]
            j <- 3
            for (i in 1:length(tau_ls)) {
                if (sum(tau_ls[[i]]) != 0 & var(tau_ls[[i]]) != 0 & !any(is.na(tau_ls[[i]])) ) {
                    tau_stack <- tau_stack + stack[j] * tau_ls[[i]]
                    j <- j + 1
                }
            }

        } else {
            print("coefs")
            print(stack)

            tau_stack <- stack[2]
            j <- 3
            for (i in 1:length(tau_ls)) {
                if (sum(tau_ls[[i]]) != 0 & var(tau_ls[[i]]) != 0 & !any(is.na(tau_ls[[i]])) ){
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
        stack_calib <- tau_calibration(tau.hat = tau_stack, Y = ho_Y, m.hat = Y_hat, W = ho_W, W.hat = 0.5)
        write.csv(stack_calib, file = paste0("./res/calib/", trial, "_Rstack_calib_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"), row.names = FALSE)

        write.csv(stack, file = paste0("./res/coeff/", trial, "_Rstack_coeffs_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"), row.names = FALSE)

        cf_tau[ho_id,] <- tau_stack # save out of bag tau estimation to dataframe
        # reformat comparison table
        method_list <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XRF", "XBART", "RBoosting", "RLasso", "RSlearner", "stack")
        rownames(compare_rloss) <- method_list
        colnames(compare_rloss) <- c("b", "c", "alpha", "mse", "mse_debiased")

        # Standardize R-loss based on outcome scale
        compare_rloss$maxmin_norm_mse <- compare_rloss$mse / (max(Y) - min(Y))
        compare_rloss$sd_norm_mse <- compare_rloss$mse / sd(Y)

        write.csv(compare_rloss, paste0("./res/", trial, "_model_sel_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"), row.names = TRUE)
    } # end of Q-fold loop
    return (cf_tau)
}

