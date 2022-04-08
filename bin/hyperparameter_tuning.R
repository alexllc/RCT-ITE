# Script containing functions that determines the hyperparameters resulting in the fewest cross-validated errors

#
# Causal Boosting parameters
#


find_cb_param <- function(x = NULL, y = NULL, tx = NULL, num_search_rounds = 5) {

    cv_performance <- data.frame()
    for(iter in 1:num_search_rounds) {

        print(paste0("Number of search round: ", iter))
        gen_param <-  list(num.trees = sample(c(100, 200, 500, 1000, 5000), 1), 
                    maxleaves = sample(c(4, 5, 10), 1), 
                    eps = sample(c(5e-6, 5e-5, 1e-4, 5e-3, 1e-2), 1), 
                    splitSpread = sample(c(0.1, 0.2, 0.3, 0.5), 1),
                    nfolds = sample(c(5, 4, 3, 2), 1)
                )
        cb <- do.call(cv.causalBoosting, append(list(x = x, y = y, tx = tx), gen_param))
        result <- c(cb$num.trees.min.effect, mean(cb$cvm.effect), gen_param)
        print("Num of min trees and avg cvm effect:")
        print(result)
        cv_performance <- rbind(cv_performance, unlist(result))
    }
    colnames(cv_performance) <- c("num_trees_min_effect", "mean_cvm_effect", "max_num_trees", "max_leaves", "eps", "split_spread", "nfolds")

    return(cv_performance)
}

#
# CausalMARS parameters
#

find_cm_param <- function(x = NULL, y = NULL, tx = NULL, validation_fold = 4, num_search_rounds = 50, verbose = FALSE) {

    cv_performance <- data.frame()

    for (iter in 1:num_search_rounds) {
        print(paste0("Searh round number: ", iter))
        gen_param <- list(maxterms = sample(c(11, min(201, max(20, 2 * ncol(x))) + 1, 3, 5, 7, 15, 21, 51, 101), 1),
                            eps = sample(c(5e-6, 5e-5, 1e-4, 5e-3, 1e-2, 0.1, 0.2, 0.5, 1), 1),
                            degree = floor(sample(c(0.1, 0.2, 0.5, 0.8, 1, 1.2, 1.5), 1) * ncol(x))
        )

        # randomly split some data for cross validation
        valID <- sample(1:dim(x)[1], floor(dim(x)[1] / validation_fold))
        cm_train_x <- x[!(1:dim(x)[1] %in% valID), ]
        cm_train_y <- y[!(1:dim(x)[1] %in% valID) ]
        cm_train_tx <- tx[!(1:dim(x)[1] %in% valID) ]

        fit <- try(cm <- do.call(causalMARS, append(list(x = cm_train_x, y = cm_train_y, tx = cm_train_tx, backstep = TRUE, x.val = x[valID,], y.val = y[valID], tx.val = tx[valID]), gen_param)))

        if (class(fit) == "try-error") {
            if (verbose) {
                print("These parameters failed: ")
                print(gen_param)
            }
           next

        } else {
            if (verbose) {
                print("These parameters worked: ")
                print(gen_param)
                print(str(fit))
            }
            cv_performance <- rbind(cv_performance, unlist(append(gen_param, sum(fit$rsstesthat))))
        }
    }
    colnames(cv_performance) = c("maxterms", "eps", "degree", "sum_rsstesthat")
    cv_performance <- cv_performance[order(cv_performance$sum_rsstesthat, decreasing = FALSE),]
    return(cv_performance) # return params that lead to smallest vaildation RSS
}

#
# PTO forest parameters
#

find_ptof_param <- function(x = NULL, y = NULL, tx = NULL, validation_fold = 4, num_search_rounds = 50) {
    
    ptof_cv_performance <- data.frame()
    for (iter in 1:num_search_rounds) {
        print(paste0("Iteration: ", iter))
        gen_param <- list(num.trees = sample(c(500, 1000, 5000, 10000, 50000, 100000), 1), 
                        mtry = sample(c(ncol(x), ceiling(min(ncol(x), sqrt(ncol(x)) + 20))), 1), 
                        min.node.size = sample(c(floor(max(25, nrow(x) / 40) * (validation_fold - 1 / validation_fold)), 25, 10, 5, 4), 1)
                        )
        ptof_error <- do.call(cv.PTOforest, append(list(x = x, y = y, tx = tx, validation_fold = validation_fold), gen_param))
        ptof_cv_performance <- rbind(ptof_cv_performance, c(unlist(gen_param), ptof_error))
    }
    colnames(ptof_cv_performance) <- c("num_trees", "mtry", "min_node_size", "mean_Z_cv_error")
    ptof_cv_performance <- ptof_cv_performance[order(ptof_cv_performance$mean_Z_cv_error),]
    
    return(ptof_cv_performance)
}