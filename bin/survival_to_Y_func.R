library(KMsurv)
library(imputeYn)
library(pseudo)

#' All in one function that accepts censoring status and time tuples and produce the estimated Y using different right-censored data handling method.
#'
#' 
#' @param event binary censoring status
#' @param time time to event
#' @param X covariate matrix in dimension n * p, optional if method of choice is "datta"
#' @param method method used to handle this right-censored data. The available methods are "datta", "datta.impYn", "pseudo", "all". Default is "datta".
#' @param impYn.method method to be used in impYn. Not required if "datta.impYn" is not chosen as the method. Default is "condMean".
#' @return A vector if a single method was chosen or a data frame with columns of output from all methods. The returned data is in the original order of your input.
#' @example surv_to_Y(time = os_years, event = censoring_status, X = covar_matrix, method = "all"))
#' @example surv_to_Y(time = os_years, event = censoring_status, X = covar_matrix, method = "datta.impYn", "condMedian"))
#' @references Datta, Somnath. “Estimating the Mean Life Time Using Right Censored Data.” Statistical Methodology 2, no. 1 (March 1, 2005): 65–69. https://doi.org/10.1016/j.stamet.2004.11.003.
#' @references Khan, Md Hasinur Rahaman, and J. Ewart H. Shaw. “On Dealing with Censored Largest Observations under Weighted Least Squares.” Journal of Statistical Computation and Simulation 86, no. 18 (December 11, 2016): 3758–76. https://doi.org/10.1080/00949655.2016.1185794.
#' @references Andersen, Per Kragh, and Maja Pohar Perme. “Pseudo-Observations in Survival Analysis.” Statistical Methods in Medical Research 19, no. 1 (February 1, 2010): 71–99. https://doi.org/10.1177/0962280209105020.


surv_to_Y <- function(event = NULL, time = NULL, X = NULL, method = "datta", impYn.method = "condMean") {

    # Kai's implementation of Datta's imputation method

    impute.survival <- function(surv.time, censor, cluters = NULL){
        # note that cluters should be a factor of length of surv.time
        if(is.null(cluters)){
            max.censored <- max(surv.time[censor == 0])
            censor[surv.time == max.censored] <- 1
            imputed.surv.times <- compute.surv.time(surv.time, censor)
        }else{
            imputed.surv.times <- rep(0, length(cluters))
            cluter_id <- levels(cluters)

            for(type in cluter_id){
                max.censored <- max(surv.time[censor == 0 & cluters == type])
                censor[surv.time == max.censored & cluters == type] <- 1
                
                idx <- (cluters == type)
                idx.surv.time <- compute.surv.time(surv.time[idx], censor[idx])
                imputed.surv.times[idx] <- idx.surv.time
            }
            return(imputed.surv.times)
        }
    }

    compute.surv.time <- function(surv.time, censor){
        # @surv.time: a vector of times, if the subject is alive, then it's NA; else, it's a number, indicated survival time
        # @censor: The status indicator, normally 0 = alive, 1 = dead, required by "Surv" function
        # need to change censoring status with the longest censoring to non-censor and assume death

        surv.obj <- Surv(surv.time, censor)
        fit.obj <- survfit(surv.obj ~ 1)
        step.surv <- -diff(c(1, fit.obj$surv))
        log.step <- log(fit.obj$time) * step.surv

        # impute the survival time of censored obs.(censor == 0), unknown to us.
        log.times <- log(surv.time)

        surv.times <- sapply(1:length(surv.time), function(i){
            ind <- which(fit.obj$time == surv.time[i])
            
            # It may be caused by a bug of survival package, sometimes we cannot find specific survival time 
            # from "time" list of fitted object by Surv function.
            if (length(ind) == 0){
                return(NA)
            }
            
            if(censor[i] == 0){
                return(sum(log.step[fit.obj$time > fit.obj$time[ind]]) / fit.obj$surv[ind])
            }else{
                return(log.times[i])
            }
        })
        return(surv.times)

    }

    if (method == "datta"){

        return(exp(impute.survival(surv.time = time, censor = event)))

    } else if(method == "datta.impYn") {
        
        if(is.null(X))
            stop("You must supply a covariate matrix to use Datta's method with imputeYn or pseudo-observation!")
        
        Yn_plus <- imputeYn(X = as.matrix(X), Y = time, delta = event, method = impYn.method)
        time[which(time == max(time))] <- Yn_plus$Yn

        return(exp(impute.survival(surv.time = time, censor = event)))
    } else if(method == "pseudo") {
        
        if(is.null(X))
            stop("You must supply a covariate matrix to use pseudo-observations!")
        
        return(pseudomean(time = time, event = event, tmax = floor(max(time))))
    } else if (method == "all") {

        all_methods <- c("datta", "datta.impYn", "pseudo") 
        output <- data.frame(id = 1:length(time))
        for (meth in methods) {
            estimate <- surv_to_Y(time = time, event = event, X = X, method = meth)
            output <- cbind(output, estimate)
            colnames(output)[which(methods %in% meth) + 1] <- meth
        }
        output$id <- NULL
        return(output)

    } else {
        stop("Invalid method! The available methods are: datta, datta.impYn, pseudo, all.")
    }
}