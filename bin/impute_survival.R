#' Function to perform 3 types of survival time imputation for right censored outcomes
#' 
impute_survival <- function(T = NULL, C = NULL, X = NULL) {
    # 1. Efron's tail correction
    source("./bin/datta_surv_imp.R")
    efron_tail <- exp(impute.survival(surv.time = T, censor = C)) # no need to convert to log scale

    # 2. Imputed tail correction (Khan & Shaw)
    # sort matrix before imputing 
    if (C[which.max(T)] == 0) {    
        chron_os <- order(T)
        imputeYn <- imputeYn(X = as.matrix(X[chron_os,]), Y = T[chron_os], delta = C[chron_os], method = "condMean") # don't impute in the log scale or else it can't handle 0s
        imputedYn <- imputeYn$Yn

        # use imputed Yn to perform Efron's tail correction again
        replace_Yn_os <- T
        replace_Yn_os[order(T, decreasing = T)[1]] <- imputedYn
        efron_tail_Yn <- exp(impute.survival(surv.time = replace_Yn_os, censor = C))
    } else {
        message("Largest OS is not censored.")
        efron_tail_Yn <- NA
    }
    # 3. Pseudo-observation
    ps_mean <- pseudomean(time = T, event = C, tmax = floor(max(T))) # negative values

    return(list(efron_tail, efron_tail_Yn, ps_mean))
}