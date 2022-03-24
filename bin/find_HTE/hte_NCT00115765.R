source("./bin/load_lib.R")
source("./bin/hyperparameter_tuning.R")
source("./bin/mod_PTOforest.R")
source("./bin/mod_causalboost.R")


trial <- "NCT00115765"

source(paste0("./bin/load_CRC_RCT/load_", trial, ".R"))
X <- as.matrix(get(trial)[[1]])
Y <- get(trial)[[2]][[1]]
W <- get(trial)[[3]]

# R-LEARNER BOOSTING has lowest MSE

Y.boost <- cvboost(X, Y, objective = "reg:squarederror", nthread = 40) # non-binary outcome
Y.hat.boost <- predict(Y.boost, newx = X)

Y.lasso <- cv.glmnet(X, Y, keep = TRUE, family = "gaussian")
Y.hat.lasso <- predict(Y.lasso, newx = X)[,1]

# which method has a smaller CV error?
print(round(c(RMSE(Y.hat.boost, Y), RMSE(Y.hat.lasso, Y)), 4))

if( RMSE(Y.hat.boost, Y) < RMSE(Y.hat.lasso, Y) ) {
    Y.hat <- Y.hat.boost
} else {
    Y.hat <- Y.hat.lasso
}

rb <- rboost(X, W, Y, p_hat = 0.5, m_hat = Y.hat, nthread = 40)
rb_tau <- predict(rb, X)


# Find if HTE is significant

test_hte <- data.frame(Y  = Y, Z = W)

hte_no_adj <- detect_idiosyncratic(formula = Y ~ Z, data = test_hte, plugin = TRUE, test.stat = "SKS.stat",  tau.hat = rb_tau, n.cores = 40)

# FRT Plug-in Test for Treatment Effect Heterogeneity 
#    Statistic P-Value (Sweep) P-Value (Plug-In)
# 1 0.06195539          0.4481            0.4441