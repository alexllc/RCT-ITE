require(grf)
require(RNOmni)
require(doParallel)
require(weights) 
# require(cwhmisc) 
require(spatstat)
require(sjstats)
require(robsurvey) 

weighted.var <- function(x, w = NULL, loc = 'mean', na.rm = FALSE) {
    # sample mean is used by default in calculation of weighted variance.
    if(!(loc %in% c('mean', 'median'))){
        stop('Error: Wrong value for loc argument!')
    }
    
    if (na.rm) {
        na <- is.na(x) | is.na(w)
        x <- x[!na]
        w <- w[!na]
    }

    if(loc == 'mean'){
        v <- sum(w * (x - weighted.mean(x, w)) ^ 2) / (sum(w))
    }else{
        v <- sum(w * (x - weighted.median(x, w)) ^ 2) / (sum(w))
    }
    return(v)
}

wtd_kurtosis <- function(x, w, loc = 'mean'){
    # n_obs <- length(w)
    if(loc == 'mean'){
        wtd_mean <- weighted.mean(x, w)
    }else if (loc == 'median'){
        wtd_mean <- weighted.median(x, w)
    }else{
        stop('Wrong value for loc argument!')
    }

    wtd_var <- weighted.var(x, w, loc = loc)
    # sample_ratio <- n_obs * (n_obs + 1)/ ((n_obs - 1) * (n_obs - 2) * (n_obs - 3))
    return((sum(w *(x - wtd_mean)^4) / sum(w)) / wtd_var^2 )
}

variance_ratio_test <- function(Y, Z, tx = 1, ctl = 0, PS = 1/2, rank_transform = T, verbose = F){
    if (rank_transform) { 
       Y <- RankNorm(Y)
    }
    Y1 <- Y[Z == tx]
    Y0 <- Y[Z == ctl]
    Pr.Z  <- mean(Z)
    wt <- ifelse(Z == tx, Pr.Z/PS, (1-Pr.Z)/(1-PS))

    N1 <- sum(wt[Z == tx])
    N0 <- sum(wt[Z == ctl])
    
    var1 <- weighted.var(Y1, wt[Z == tx], loc = 'median') 
    var0 <- weighted.var(Y0, wt[Z == ctl], loc = 'median') 

    log.varR <- log(var1/var0)

    asy.se <- sqrt(abs((wtd_kurtosis(Y1, wt[Z == tx]) - 1)/N1 + (wtd_kurtosis(Y0, wt[Z == ctl]) - 1)/N0))
    pvalue <- as.numeric((1 - pnorm(abs(log.varR), 0, asy.se)))*2

    res <- data.frame(pvalue = pvalue, var1 = var1, var0 = var0)
    res$ratio = res$var1/res$var0
    res$log.ratio = log(res$ratio)
    res$asy.se = asy.se
    res$z = res$log.ratio / res$asy.se
    if(verbose){
        cat('kurtosis pvalue:', pvalue, fill = T)
    }
    return(c(kurtosis_stat = res$z, kurtosis_pval = pvalue))
}

weighted_Levene_test <- function(Y, Z, tx = 1, ctl = 0, PS = 1/2, 
                                 rank_transform = T,
                                 location = "median", 
                                 trim.alpha = 0.1,
                                 verbose = TRUE)  {
    # @param Y The outcome for the test.
    # @param Z Treatment assignment variable.
    # @param tx treatment group indicator.
    # @param ctl Control group indicator.
    # @param PS Proposity scores, should be of the same length as Outcome.
    # @param rank_transform Indicate whether rank transformation of outcome is applied.
    # @param location The method for computation of mean effect, "median" is used by default.
    # @param trim.alpha The proportion of extremed observations trimed away.
    # @param verbose The switch for output.

    if(!(location %in% c('mean', 'median', 'trim.mean'))){
        stop('Error: Wrong value for location argument!')
    }

    if (rank_transform) { 
       Y = RankNorm(Y)
    }
   
    Pr.Z  <- mean(Z)
    wt <- ifelse(Z == tx, Pr.Z/PS, (1-Pr.Z)/(1-PS))
    
    if(any(wt) < 0){
        return(rep(NA,6))
    }
    Y1 <- Y[Z == tx]
    Y0 <- Y[Z == ctl]
    
    wt1 <- wt[Z == tx]
    wt0 <- wt[Z == ctl]

    N1 <- sum(wt[Z == tx])
    N0 <- sum(wt[Z == ctl])
    
    if(location == "mean"){
       mu1 = weighted.mean(Y1, wt[Z == tx]) 
       mu0 = weighted.mean(Y0, wt[Z == ctl]) 
    }else if(location == "median"){
       mu1 <- weighted.median(Y1, wt[Z == tx]) 
       mu0 <- weighted.median(Y0, wt[Z == ctl]) 
    }else {
       mu1 <- weighted_mean_trimmed(Y1, wt[Z == tx], LB = trim.alpha,  na.rm = FALSE) 
       mu0 <- weighted_mean_trimmed(Y0, wt[Z == ctl], LB = trim.alpha,  na.rm = FALSE)
    }
    
    Y1_star <- abs(Y1 - mu1) 
    Y0_star <- abs(Y0 - mu0)
    
    t_statistics_obj <- wtd.t.test(Y1_star, Y0_star, weight=wt1, weighty=wt0, samedata = F,
                                   alternative ="two.tailed") 
    
    t_stat <- t_statistics_obj$coefficients[1]
    t_pval <- t_statistics_obj$coefficients[3]
  
    # or use https://lindeloev.github.io/tests-as-linear/ section 5.2
    Y_star <- sapply(seq(length(Z)), function(i) abs(ifelse(Z[i] == tx, Y[i] - mu1, Y[i] - mu0)))

    data <- data.frame(Y_star, Z, wt)
    wtd_mwu_obj <- weighted_mannwhitney(data, Y_star, Z, wt)
    wtd_mwu_statistic <- wtd_mwu_obj$statistic 
    wtd_mwu_pvalue <- wtd_mwu_obj$p.value
    
    lm_obj <- lm(Y_star ~ Z, weights = wt)  
    lm_t_stat = summary(lm_obj)$coef[2,3]
    lm_pval = summary(lm_obj)$coef[2,4]

    if(verbose){
        cat('pvalue of t.test:', t_pval, fill = T)
        cat('pvalue of wtd_mwu:', wtd_mwu_pvalue, fill = T) # Mann-Whitney U test 
    }
        
    return(c(mwu_stat = wtd_mwu_statistic, mwu_pval = wtd_mwu_pvalue,  
             t_stat = t_stat, t_pval = t_pval,
	     lm_stat = lm_t_stat, lm_pval = lm_pval))
}
    
