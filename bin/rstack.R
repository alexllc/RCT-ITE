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

perform_rstack <- function(Y = NULL, X = NULL, W = NULL, trial_name = NULL, outcome = NULL, imp_type_Y = NULL, prob = 0.5, Q = 4, tuned_cb_param = FALSE, tuned_cm_param = TRUE, tune_ptof_param = TRUE, perform_xb = FALSE) {

    available_sample <- 1:dim(X)[1]
    cf_tau <- data.frame(ID = available_sample)
    rownames(cf_tau) <- cf_tau$ID

    for (iter in 1:Q) {

        tau_ls <- list()

        if (iter != Q) {
            ho_id <- sample(available_sample, floor(dim(X)[1] / Q))
            available_sample <- available_sample[!(available_sample %in% ho_id)]
        } else {
            ho_id <- available_sample
        }
        train_id <- seq(1:dim(X)[1])[!(seq(1:dim(X)[1]) %in% ho_id)]

        train_X <- X[train_id,]
        train_Y <- Y[train_id]
        train_W <- W[train_id]

        ho_X <- X[ho_id,]
        ho_Y <- Y[ho_id]
        ho_W <- W[ho_id]

        #
        # R-stack
        #
        mhat <- get_mhat(X_train = train_X, W_train = train_W, Y_train = train_Y, X_ho = ho_X, Y_ho = ho_Y)
        
        #
        # Causal forest
        #
        message(paste0(rep("=", 80)))
        message(paste0("Building Causal Forest."))
        message(paste0(rep("=", 80)))

        cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
        cf_tau_pred <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)
        while (var(cf_tau_pred$predictions) == 0) { # tuned forests sometimes give equal predictions, if that's the case, rebuild forest
            cf <- causal_forest(train_X, train_Y, train_W, num.trees = 100000, tune.parameters = "all")
            cf_tau_pred <- predict(cf, newdata = ho_X, estimate.variance =  TRUE)
        }


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
        

        #
        # CausalMARS
        #
        message(paste0(rep("=", 80)))
        message(paste0("Building Performing causal MARS."))
        message(paste0(rep("=", 80)))
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
            write.csv(cm_cv_performance, file = paste0(paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_cm_params_used.csv")))

        } else {
            cm <- causalMARS(x = train_X, tx = train_W, y = train_Y)
        }


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
            if (any(is.na(ptof_tau_pred))) {
                message("Out of bag NA predictions generated, will be refitting the model.")
                rebuild_ptof <- any(is.na(ptof_tau_pred))
            } else {
                rebuild_ptof <- FALSE
            }
        }
        write.csv(ptof_param_search, file = paste0(paste0("./res/params/", trial, "_", imp_type_Y, "_", outcome, "_ptof_params_used.csv")))

        if (class(ptof) == "try-error" | class(ptof_tau_pred) == "try-error") {
            print("Cannot build this model.")
            ptof_tau_pred <- rep(0, dim(ho_X)[1])
        }


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
        } else {
            xb_tau_pred <- rep(0, dim(ho_X)[1])
        }
        
        #
        # X-learner of random forests
        #

        x_rf <- X_RF(feat = train_X, tr = train_W, yobs = train_Y)
        x_rf_tau <- EstimateCate(x_rf, ho_X)

        #
        # X-learner of BART
        #

        x_bart <- X_BART(feat = train_X, tr = train_W, yobs = train_Y)
        x_bart_tau <- EstimateCate(x_bart, ho_X)


        ##################################################
        # R-learner
        ##################################################
        message(paste0(rep("=", 80)))
        message(paste0("Building R-learner."))
        message(paste0(rep("=", 80)))

        #
        # R-learner of boosting
        #

        rb <- try(rboost(train_X, train_W, train_Y, p_hat = prob, m_hat = mhat, nthread = 40))
        while(class(rb) == "try-error") {
            rb <- try(rboost(train_X, train_W, train_Y, p_hat = prob, m_hat = mhat, nthread = 40))
        }
        rb_tau <- predict(rb, ho_X)

        #
        # R-learner of LASSO
        #

        rl <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = mhat))
        while(class(rl) == "try-error") {
            rl <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = mhat))
        }
        rl_tau <- predict(rl, ho_X)[,1]


        #
        # RS-learner of LASSO
        #

        rl_RS <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = mhat, rs = TRUE)) # returned same predictions for all patients
        while(class(rl_RS) == "try-error") {
            rl_RS <- try(rlasso(train_X, train_W, train_Y, p_hat = prob, m_hat = mhat, rs = TRUE))
        }
        rl_RS_tau <- predict(rl_RS, ho_X)[,1]


        tau_ls <- append(tau_ls, list(cf_tau_pred$predictions, cb_tau_pred, tau_cm_pred, ptof_tau_pred, xb_tau_pred, x_rf_tau, x_bart_tau, rb_tau, rl_tau, rl_RS_tau))

        names(tau_ls) <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XRF", "XBART", "RBoosting", "RLasso", "RSlearner")


        ##################################################
        # R-stack
        ##################################################

        RESP <- ho_Y - mhat
        R_mat <- cbind(1, ho_W - rep(prob, length(ho_W)))
        tau_counter <- 0
        for(tau_pred in tau_ls) {
            # nnls will fail due to a non-positive definite matrix if you have a constant tau vector
            # sometimes tree based tau estimators will render NAs in out of bag tau predictions, which will be excluded in the stack
            if (sum(tau_pred, na.rm = TRUE) != 0 & var(tau_pred) != 0 & !any(is.na(tau_pred))){ 
            
                R_mat <- cbind(R_mat, (ho_W - rep(prob, length(ho_W))) * tau_pred)
                tau_counter <- tau_counter + 1
            }
        }

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
            mean((ho_Y - mhat - tau_stack * (ho_W - rep(prob, length(ho_W))))^2)

            tau_stack_rloss <- try(rloss(tau.pred = tau_stack, Y = ho_Y, W = ho_W, Y.hat = mhat, prob = prob))

        }

        cf_tau[ho_id,] <- tau_stack # save out of bag tau estimation to dataframe
    }
    return (cf_tau)
}

