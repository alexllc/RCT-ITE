clin_df <- read.csv("./dat/NCT00364013_adls.csv")
biom <- read.csv("./dat/NCT00364013_biom.csv")

# Subjects with at least one biomarker entry is excluded from model building
tested_subj <- Reduce(union, lapply(biom[,c(2:dim(biom)[2])], function(X) which(is.na(X))))
test_clin_df <- clin_df[clin_df$SUBJID %in% tested_subj,]

dim(test_clin_df)

test_clin_df$PFS_mo <- test_clin_df$PFSDYCR / 30.4167
test_clin_df$OS_mo <- test_clin_df$DTHDY / 30.4167

blls <- read.csv("./dat/NCT00364013_blls.csv")
tlls <- read.csv("./dat/NCT00364013_tlls.csv")
rsp <- read.csv("./dat/NCT00364013_rsp.csv")


bigX <- left_join(test_clin_df, biom, by = c("SUBJID"))
bigX <- left_join(bigX, blls, by = c("SUBJID"))
bigX <- left_join(bigX, tlls, by = c("SUBJID"))
bigX <- left_join(bigX, rsp, by = c("SUBJID"))


X <- select(test_clin_df, all_of(c("LIVERMET", "DIAGMONS", "AGE",  "SEX", "B_WEIGHT", "B_HEIGHT", "RACE",  "B_ECOG", "HISSUBTY", "B_METANM", "DIAGTYPE")))
W <- test_clin_df$ATRT

########################################################################
#### Step 2: address outcome censoring
########################################################################

# 1. Efron's tail correction
source("./audit/original_surv.r")
efron_tail <- exp(impute.survival(surv.time = test_clin_df$PFS_mo, censor = test_clin_df$PFSCR)) # no need to convert to log scale

# 2. Imputed tail correction (Khan & Shaw)
# sort matrix before imputing 
chron_os <- order(test_clin_df$PFS_mo)
imputeYn <- imputeYn(X = as.matrix(X[chron_os,]), Y = log(test_clin_df$PFS_mo[chron_os]), delta = test_clin_df$PFSCR[chron_os], method = "condMean")

# largest observation is not censored, no need to impute Yn

# use imputed Yn to perform Efron's tail correction again
efron_tail <- exp(impute.survival(surv.time = test_clin_df$PFS_mo, censor = test_clin_df$PFSCR))

# 3. Pseudo-observation
ps_mean <- pseudomean(time = test_clin_df$PFS_mo, event = test_clin_df$PFSCR, tmax = floor(max(test_clin_df$PFS_mo))) # negative values



# change names from source data
X_test <- as.matrix(X)
Y_test <- efron_tail
W_test <- W

prob <- 0.5 # randomized experiment

compare_rloss <- data.frame(b = numeric(), c = numeric(), alpha = numeric(), mse = numeric(), debiased_mse = numeric())


# Since CausalBoost is not implemented w.r.t. the causal forest in Imbens and Athey, there is not debiased.error to take out of context

# let's first fit cvlasso or cv boosting to estimate nuisance component marginal response model (m.hat)

Y.boost <- cvboost(X_test, Y_test, objective = "reg:squarederror", nthread = 4) # non-binary outcome
Y.hat.boost <- predict(Y.boost)

Y.lasso <- cv.glmnet(X_test, Y_test, keep = TRUE, family = "gaussian")
Y.hat.lasso <- Y.lasso$fit.preval[,!is.na(colSums(Y.lasso$fit.preval))]
Y.hat.lasso <- Y.hat.lasso[, Y.lasso$lambda == Y.lasso$lambda.min]

print(round(c(RMSE(Y.hat.boost, Y), RMSE(Y.hat.lasso, Y)), 4))
# boost had a smaller CV error
Y.hat <- Y.hat.boost

ptof_test <- PTOforest(x = X_test, tx = W_test, y = Y_test, num.trees = 10000)
tau_ptof_pred_test <- predict(ptof, X_test)
ptof_rloss_pred_test <- rloss(tau.pred = tau_ptof_pred_test, Y = Y_test, W = W_test, Y.hat = Y.hat, prob = prob) # constant?

# adding baseline additional covar
bigX_test <- left_join(test_clin_df, blls, by = c("SUBJID"))


# convert categorical to numeric
num_bigX_test <- cat2num(bigX_test)
num_bigX <- num_bigX_test[[1]]
bigX_ddt <- num_bigX_test[[2]]

#
# impute missing data if necessary
#
if (any(is.na(num_bigX))) {
    pre_proc <- scale(num_bigX, center = TRUE, scale = TRUE)
    imp_bigX <- impute.knn(pre_proc, k=10)
    message(c("Seed used: ", imp_bigX$rng.seed))
    imp_bigX <- t(apply(imp_bigX$data, 1, function(r) r *attr(imp_bigX$data,'scaled:scale') + attr(imp_bigX$data, 'scaled:center')))
}
message("Is there still missing data? ", any(is.na(imp_bigX)))
imp_bigX <- as.data.frame(imp_bigX)
class(imp_bigX$PFSCR) <- "integer" # after imputation returned floating point for 0s

bl_bigX <- select(test_clin_df, -all_of(c("SUBJID", "ATRT", "PRSURG", "DTHDY", "DTH", "PFSDYCR", "PFSCR", "PFS_mo", "OS_mo")))
imp_bigX_Y <- exp(impute.survival(surv.time = bl_bigX$PFS_mo, censor = bl_bigX$PFSCR)) # no need to convert to log scale

ptof_test_bl <- PTOforest(x = bl_bigX, tx = imp_bigX$ATRT, y =imp_bigX, num.trees = 10000)
tau_ptof_pred_test <- predict(ptof, imp_bigX)
ptof_rloss_pred_test <- rloss(tau.pred = tau_ptof_pred_test, Y = Y_test, W = W_test, Y.hat = Y.hat, prob = prob) # constant?


plot_df <- data.frame(SUBJID = test_clin_df$SUBJID, tau = tau_ptof_pred_test, sign = ifelse(tau_ptof_pred_test > 0, "positive", "negative"))
plot_df <- left_join(plot_df, biom, by = c("SUBJID"))

ggbarplot(plot_df, x = "SUBJID", y = "tau",
          fill = "sign",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in ascending order
          sort.by.groups = TRUE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "treatment effect (panitumumab + FOLFOX) - FOLFOX alone",
          xlab = FALSE,
          legend.title = "TE"
          )

rl_test <- rlasso(X_test, W_test, Y_test, p_hat = prob, m_hat = Y.hat)


library(survminer)
library(survival)
fit <- survfit(Surv(PFS_mo, PFSCR) ~ LIVERMET,
               data = test_clin_df)
# Visualize with survminer
ggsurvplot(
   fit,                     # survfit object with calculated statistics.
   data = test_clin_df,  # data used to fit survival curves. 
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   conf.int = TRUE,         # show confidence intervals for 
                            # point estimaes of survival curves.
                            # survival estimates.
   ggtheme = theme_minimal(), # customize plot and risk table with a theme.
 risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
)

extreme_sub <- test_clin_df$SUBJID[order(plot_df$tau)[1:4]]