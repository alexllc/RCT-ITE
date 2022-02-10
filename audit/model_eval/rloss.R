source("./audit/model_eval/load_lib.R")
source("./audit/model_eval/load_data.R")

rloss <- function(tau.pred = NULL, Y = NULL, W = NULL, Y.hat = NULL, prob = NULL) {
    RESP = Y - Y.hat
    R.mat = cbind(1, W - prob,
                    (W - prob) * tau.pred)

    learn_coeff = nnls(R.mat, RESP, constrained = c(FALSE, FALSE, TRUE))

    print("coefs")
    print(learn_coeff)
    debiased_tau <- learn_coeff[2] + learn_coeff[3] * tau.pred
    mse_debiased <- mean((Y - Y.hat - debiased_tau * (W - prob))^2)
    mse <- mean((Y - Y.hat - tau.pred * (W - prob))^2)

    return(list(nnls_coeff = learn_coeff, debiased_tau = debiased_tau, mse_debiased = mse_debiased, mse = mse))
}

##################################################
# Fit various ITE mdoels and calculate R-loss
##################################################

# change names from source data
X <- as.matrix(X)
Y <- efron_tail
W <- W

prob <- 0.5 # randomized experiment

compare_rloss <- data.frame(b = numeric(), c = numeric(), alpha = numeric(), mse = numeric(), debiased_mse = numeric())


# Since CausalBoost is not implemented w.r.t. the causal forest in Imbens and Athey, there is not debiased.error to take out of context

# let's first fit cvlasso or cv boosting to estimate nuisance component marginal response model (m.hat)

Y.boost <- cvboost(X, Y, objective = "reg:squarederror", nthread = 4) # non-binary outcome
Y.hat.boost <- predict(Y.boost)

Y.lasso <- cv.glmnet(X, Y, keep = TRUE, family = "gaussian")
Y.hat.lasso <- Y.lasso$fit.preval[,!is.na(colSums(Y.lasso$fit.preval))]
Y.hat.lasso <- Y.hat.lasso[, Y.lasso$lambda == Y.lasso$lambda.min]

print(round(c(RMSE(Y.hat.boost, Y), RMSE(Y.hat.lasso, Y)), 4))
# lasso had a smaller CV error
Y.hat <- Y.hat.lasso

#
# Causal forest
#

# R-loss already implemented by the package and is contained in the debiased.error column in predict.causal_forest(). You can obtain MSE by mean(cf$debiased.error), the values were already squared.
# since causal trees are honest already, there's no need to futher split into k-fold estimation

cf <- causal_forest(X, Y, W)
cf_tau_pred <- predict(cf, estimate.variance = TRUE)
print(mean(cf_tau_pred$debiased.error))
# equivalent to
cf_mse <- mean((Y - Y.hat - cf_tau_pred$predictions * (W - prob))^2)
cf_rloss_res <- rloss(tau.pred = cf_tau_pred$predictions, Y = Y, W = W, Y.hat = Y.hat, prob = prob)

compare_rloss <- rbind(compare_rloss, c(cf_rloss_res$nnls_coeff, cf_rloss_res$mse, cf_rloss_res$mse_debiased))


##################################################
# Powers: CausalBoost, CausalMARS, PTOForest
##################################################

# now we estimate tau from CausalBoosting

cb <- causalBoosting(x = X, tx = W, y = Y) # needs hyperparameter tuning
min_y_err <- cb$err.y[which(cb$err.y == min(cb$err.y))]
tau_cb_pred <- cb$tauhat[,which(cb$err.y == min(cb$err.y))] # choosing tau that has the minimum y.err

cb_rloss_res <- rloss(tau.pred = tau_cb_pred, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(cb_rloss_res$nnls_coeff, cb_rloss_res$mse, cb_rloss_res$mse_debiased))


#
# CausalMARS
#

cm <- causalMARS(x = X, tx = W, y = Y)
tau_cm_pred <- predict(cm, X)
cm_rloss_pred <- rloss(tau.pred = tau_cm_pred, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(cm_rloss_pred$nnls_coeff, cm_rloss_pred$mse, cm_rloss_pred$mse_debiased))


#
# PTOForest
#
ptof <- PTOforest(x = X, tx = W, y = Y)
tau_ptof_pred <- predict(ptof, X)
ptof_rloss_pred <- rloss(tau.pred = tau_ptof_pred, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(ptof_rloss_pred$nnls_coeff, ptof_rloss_pred$mse, ptof_rloss_pred$mse_debiased))



##################################################
# X-learner
##################################################
xb <- xboost(x = X, w = W, y = Y)
tau_xb_pred <- predict(xb, X)

xb_rloss <- rloss(tau.pred = tau_xb_pred, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(xb_rloss$nnls_coeff, xb_rloss$mse, xb_rloss$mse_debiased))


x_rf <- X_RF(feat = X, tr = W, yobs = Y)
x_rf_tau <- EstimateCate(x_rf, X)
x_rf_rloss <- rloss(tau.pred = x_rf_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(x_rf_rloss$nnls_coeff, x_rf_rloss$mse, x_rf_rloss$mse_debiased))


x_bart <- X_BART(feat = X, tr = W, yobs = Y)
x_bart_tau <- EstimateCate(x_bart, X)
x_bart_rloss <- rloss(tau.pred = x_bart_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(x_bart_rloss$nnls_coeff, x_bart_rloss$mse, x_bart_rloss$mse_debiased))


##################################################
# R-learner
##################################################
rb <- rboost(X, W, Y, p_hat = prob, m_hat = Y.hat, nthread = 4)
rb_tau <- predict(rb)

rl <- rlasso(X, W, Y, p_hat = prob, m_hat = Y.hat)
rl_tau <- predict(rl)

print(round(c(mean((Y - Y.hat - rb_tau * (W - prob))^2), mean((Y - Y.hat - rl_tau * (W - prob))^2)), 4)) # boosting wins

rb_rloss <- rloss(tau.pred = rb_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(rb_rloss$nnls_coeff, rb_rloss$mse, rb_rloss$mse_debiased))

# RS learner
rl_RS <- rlasso(X, W, Y, p_hat = prob, m_hat = Y.hat, rs = TRUE) # returned same predictions for all patients
rl_RS_tau <- predict(rl_RS)[,1]
rs_rloss <- rloss(tau.pred = rl_RS_tau, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
# do not include

##################################################
# R-stack
##################################################

RESP = Y - Y.hat
R.mat <- cbind(1, W - prob,
                (W - prob) * cf_tau_pred$predictions,
                (W - prob) * tau_cb_pred,
                (W - prob) * tau_cm_pred,
                (W - prob) * tau_ptof_pred,
                (W - prob) * tau_xb_pred,
                (W - prob) * x_bart_tau,
                (W - prob) * x_rf_tau,
                (W - prob) * rb_tau
                )

stack = nnls(as.matrix(R.mat), RESP, constrained = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
print("coefs")
print(stack)

tau.stack = stack[2] +
    stack[3] * cf_tau_pred$predictions +
    stack[4] * tau_cb_pred +
    stack[5] * tau_cm_pred +
    stack[6] * tau_ptof_pred +
    stack[7] * tau_xb_pred
    stack[8] * x_bart_tau
    stack[9] * x_rf_tau
    stack[10] * rb_tau

mean((Y - Y.hat - tau.stack * (W - prob))^2)

tau_stack_rloss <- rloss(tau.pred = tau.stack, Y = Y, W = W, Y.hat = Y.hat, prob = prob)
compare_rloss <- rbind(compare_rloss, c(tau_stack_rloss$nnls_coeff, tau_stack_rloss$mse, tau_stack_rloss$mse_debiased))


# reformat comparison table
method_list <- c("CF", "CBoost", "CMARS", "PTOForest", "XBoosting", "XBART", "XRF", "RBoosting", "stack")
rownames(compare_rloss) <- method_list
colnames(compare_rloss) <- c("b", "c", "alpha", "mse", "mse_debiased")

write.csv(compare_rloss, "./res/NCT00364013_model_sel_res.csv", row.names = TRUE)