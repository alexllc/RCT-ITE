# Script containing functions that determines the hyperparameters resulting in the fewest cross-validated errors

#
# Causal Boosting
#
library(foreach)
library(doParallel)



find_cb_param <- function(x = NULL, y = NULL, tx = NULL, num_search_rounds = 5) {

    nCores <- detectCores() - 1
    cluster <- makeCluster(nCores)

    doParallel::registerDoParallel(cluster)

    # parallel computing version

    res <- foreach (iter = 1:num_search_rounds, .combine = rbind) %do% {

        cb_param <-  list(num.trees = sample(c(500, 1000, 10000, 50000), 1), 
                            maxleaves = sample(c(4, 5, 10), 1), 
                            eps = sample(c(5e-5, 1e-4, 5e-3, 1e-2), 1), 
                            splitSpread = sample(c(0.1, 0.2, 0.3, 0.5), 1)
                        )

        cb <- do.call(cv.causalBoosting, append(list(x = x, y = y, tx = tx), cb_param))
        result <- c(cb$num.trees.min.effect, mean(cb$cvm.effect))

        result
    }
    which(min(result[,2]))
    return(res)
    stopImplicitCluster()
}

find_cb_param <- function(x = NULL, y = NULL, tx = NULL, num_search_rounds = 5) {

    cv_performance <- list()
    for(iter in 1:num_search_rounds) {

        print(c("Number of search round: ", iter))
        gen_param <-  list(num.trees = sample(c(500, 1000, 10000, 50000), 1), 
                    maxleaves = sample(c(4, 5, 10), 1), 
                    eps = sample(c(5e-6, 5e-5, 1e-4, 5e-3, 1e-2), 1), 
                    splitSpread = sample(c(0.1, 0.2, 0.3, 0.5), 1)
                )
        cb <- do.call(cv.causalBoosting, append(list(x = x, y = y, tx = tx), gen_param))
        result <- list(cb$num.trees.min.effect, mean(cb$cvm.effect), gen_param)
        print("Num of min trees and avg cvm effect:")
        print(result)
        cv_performance <- append(gen_param, result)
}
    return(cv_performance)
}