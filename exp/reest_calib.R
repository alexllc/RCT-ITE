source("./bin/load_lib.R")
source("./bin/rstack.R")

Q <- 4

# Extract trial names from list of scripts
trial_scripts <- list.files("./bin/load_RCT")
trial_ls <- trial_scripts[grep("^load*", trial_scripts)]
trial_ls <- unlist(lapply(strsplit(trial_ls, "_|\\."), function(X) X[2]))
reest_trials <- c("NCT00003299", "NCT00041119chemo", "NCT00113763", "NCT00115765", "NCT00119613")
completed_trials <- c("NCT00041119length", "NCT00339183", "NCT00364013", "NCT00460265")

    for (trial in reest_trials) {
        trial_avg_calib_list <- list()
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
        outcome_list <- get(paste0(trial, "_outcomes"))

        for (outcome in outcome_list) {

            message(paste0("Processing outcome: ", outcome))

            Y_list <- get(paste0(outcome, "_Y_list"))

            imp_type_Y <- "efronYn"
            Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

            if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
                Y <- as.numeric(Y_list[[1]])
                imp_type_Y <- "efron"
            }

            # Load used methods
            ite_res_files <- list.files("./res/indiv_model_tau/")
            trial_res <- ite_res_files[grep(paste0("^", trial, "*"), ite_res_files)]
            tau_method_used <- unique(unlist(lapply(strsplit(trial_res, "_|\\."), function(X) X[2])))

            if (file.exists(paste0("./res/ite_tau_estimates/", trial, "_", imp_type_Y, "_", outcome, "_Rstack_tau_estimates.csv"))) {
                tau_method_used <- c(tau_method_used, "Rstack")
            }

            # Load CV folds
            ho_matrix <- read.csv(paste0("./dat/cv_ho_ids/est_tau_", trial, "_holdout_id.csv"))
            #
            # R-stack
            #
            if (!(imp_type_Y == "efronYn" | imp_type_Y == "efron")) {
                binary_Y <- TRUE
            } else {
                binary_Y <- FALSE
            }

            for (tau_meth in tau_method_used) {
                message(paste0("Processing method: ", tau_meth))
                avg_calib_list <- list()
                for (iter in 1:Q) {

                    message(paste0("Processing iteration: ", iter))
                    # Load tau_estimates
                    if (tau_meth == "Rstack") {
                        tau_pred <- read.csv(paste0("./res/ite_tau_estimates/", trial, "_", imp_type_Y, "_", outcome, "_Rstack_tau_estimates.csv"))
                    } else {
                        tau_pred <- try(read.csv(paste0("./res/indiv_model_tau/", trial, "_", tau_meth, "_", imp_type_Y, "_", outcome, "_", iter, "_tau.csv")))
                    }
                    if(class(tau_pred) == "try-error") {
                        next
                    }
                    tau_pred <- tau_pred[,1]

                    # orig_calib <- read.csv(paste0("./res/calib/", trial, "_", tau_meth, "_calib_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"))

                    # Divide datasets in to CV folds
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

                    message("ESTIMATING Y.HAT.")
                    Y_hat <- get_mhat(X_train = train_X, W_train = train_W, Y_train = train_Y, X_ho = ho_X, Y_ho = ho_Y, binary_Y = binary_Y, fast = TRUE)
                    
                    message("NEW CALIBRATION RESULTS: ")
                    # Restimate calibration
                    if (tau_meth == "Rstack") {
                        new_calib <- tau_calibration(tau.hat = tau_pred[ho_id], Y = ho_Y, m.hat = Y_hat, W = ho_W, W.hat = 0.5)
                    } else {
                        new_calib <- tau_calibration(tau.hat = tau_pred, Y = ho_Y, m.hat = Y_hat, W = ho_W, W.hat = 0.5)
                    }
                    avg_calib_list <- append(avg_calib_list, list(c(new_calib[1,], new_calib[2,])))
                } # end of iteration loop

                if(class(tau_pred) == "try-error") {
                        next
                }

                avg_calib <- do.call(rbind.data.frame, avg_calib_list)
                avg_calib$term <- NULL
                avg_calib$`term.1` <- NULL
                trial_avg_calib_list <- append(trial_avg_calib_list, list(c(outcome, tau_meth, as.list(colMeans(avg_calib)))))
            } # end of methods loop
        } # end of outcome loop

        trial_avg_calib <- do.call(rbind.data.frame, trial_avg_calib_list)
        colnames(trial_avg_calib) <- c("outcome", "tau_method", "ate.estimate", "ate.std.error", "ate.statistic", "ate.pval", "hte.estimate", "hte.std.error", "hte.statistic", "hte.pval")

        write.csv(trial_avg_calib, paste0("./res/calib/calib_avg/", trial, "_calib_res.csv"), row.names = FALSE)
    }