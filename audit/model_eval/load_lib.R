library(dplyr)

# Survival analysis
library(survival)
library(survminer)
library(imputeYn)
library(pseudo)

# Causal inference
library(causalToolbox) # X-learner
library(grf)
library(causalLearning) # Powers

# R-learner
library(rlearner)
library(KRLS2)
library(glmnet)
library(nnls)

library(caret) # cross validation

# Various causal inference framework
library(causalLearning)
library(causalToolbox)


nnls = function(M, v, constrained) {
    Dmat = t(M) %*% M
    dvec = t(M) %*% v
    Amat = matrix(0, ncol(M), sum(constrained))
    cons = which(constrained)
    for(iter in 1:length(cons)) {
        Amat[cons[iter], iter] = 1
    }
    bvec = rep(0, sum(constrained))
    soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec) 
    soln$solution
}