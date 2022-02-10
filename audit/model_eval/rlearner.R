source("./audit/model_eval/load_lib.R")
source("./audit/model_eval/load_data.R")


#############################################################################
###### RAW CODE FOLLOWS ######
#############################################################################

rlasso_fit = rlasso(as.matrix(covar), W, log(efron_tail_Yn))
rlasso_est = predict(rlasso_fit, as.matrix(covar))

rboost_fit = rboost(as.matrix(covar), W, efron_tail_Yn)
rboost_est = predict(rboost_fit, as.matrix(covar))

rkern_fit = rkern(as.matrix(covar), W, efron_tail_Yn)
rkern_est = predict(rkern_fit, as.matrix(covar))

# holdout data
testid <- sample(dim(X)[1],dim(X)[1]*0.25)
trainid <- seq(1:dim(X)[1])[!(seq(1:dim(X)[1]) %in% testid)]

X <- as.matrix(covar)[trainid,]
X.holdout <- as.matrix(covar)[testid,]

Y <- efron_tail_Yn[trainid]
Y.holdout <- efron_tail_Yn[testid]

W <- tx[trainid]
W.holdout <- tx[testid]

## Two step R learner

#
# fit propensity model
#

# do not fit propensity if randomized
W.hat <- 0.5

#
# fit marginal response model
#
Y.boost = cvboost(X, Y, objective = "reg:squarederror", nthread = 4) # non-binary outcome
Y.hat.boost = predict(Y.boost)

Y.lasso = cv.glmnet(X, Y, keep = TRUE, family = "gaussian")
Y.hat.lasso = Y.lasso$fit.preval[,!is.na(colSums(Y.lasso$fit.preval))][, Y.lasso$lambda == Y.lasso$lambda.min]

# calculate cross validated RMSE
print(round(c(RMSE(Y.hat.boost, Y), RMSE(Y.hat.lasso, Y)), 4))
Y.hat = Y.hat.boost # boosting has smaller RMSE


#
# fit R-learner given chosen nuisance components
#
tau.boost = rboost(X, W, Y, p_hat = W.hat, m_hat = Y.hat, nthread = 4)
tau.hat.boost = predict(tau.boost)


tau.lasso = rlasso(X, W, Y, p_hat = W.hat, m_hat = Y.hat)
tau.hat.lasso = predict(tau.lasso)

# who wins on CV error?
print(round(c(mean((Y - Y.hat - tau.hat.boost * (W - W.hat))^2), mean((Y - Y.hat - tau.hat.lasso * (W - W.hat))^2)), 4)) #lasso wins


# test using holdout data
tau.hat.boost.holdout = predict(tau.boost, X.holdout)
tau.hat.lasso.holdout = predict(tau.lasso, X.holdout)

Y.hat.holdout = predict(Y.boost, X.holdout)
W.hat.holdout = rep(W.hat, length(testid))

print(round(c(mean((Y.holdout - Y.hat.holdout - tau.hat.boost.holdout * (W.holdout - W.hat.holdout))^2), mean((Y.holdout - Y.hat.holdout - tau.hat.lasso.holdout * (W.holdout - W.hat.holdout))^2)), 4)) # lasso wins

