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

########################################################################
#### Step 1: Prepare data
########################################################################
df <- read.csv("/home/alex/Documents/lab/ctci/dat/PDS/Multiple_Allianc_2002_213_NCT00052910_esophageal/NCT00052910_D1_(EFFICACY).csv")
tx <- df$TREAT_ASSIGNED - 1
outcome_col <- c("death_final", "os_year", "DFS_status", "DFS_year")
outcome <- dplyr::select(df, all_of(outcome_col))
covar <- dplyr::select(df, -c("MASK_ID", "TREAT_ASSIGNED", outcome_col))

# prepare imputed covariate matrix (imported from Python)
covar <- read.csv("./audit/imputed_covar.csv")
covar$os_year <- NULL
covar$TREAT_ASSIGNED <- NULL

########################################################################
#### Step 2: address outcome censoring
########################################################################

# 1. Efron's tail correction
source("./audit/original_surv.r")
efron_tail <- exp(impute.survival(surv.time = df$os_year, censor = df$death_final)) # no need to convert to log scale

# 2. Imputed tail correction (Khan & Shaw)
imputeYn <- imputeYn(X = as.matrix(covar), Y = log(df$os_year), delta = df$death_final, method = "condMean") # imputed values not greater

# 3. Pseudo-observation
ps_mean <- abs(pseudomean(time = df$os_year, event = df$death_final, tmax = floor(max(df$os_year)))) # negative values


########################################################################
#### Step 3: Estimate ITE
########################################################################
Y_choice <- ps_mean

############################# X-learner ################################

# randomly sample half of the data without replacement as train set
set.seed(123) # set seed for reproduciblity
train_set <- sample(df$MASK_ID, size = floor(dim(df)[1]/2))
test_set <- df$MASK_ID[!df$MASK_ID %in% train_set]

xl_rf <- X_RF(feat = covar[train_set,], tr = tx[train_set], yobs = Y_choice[train_set])
xl_bart <- X_BART(feat = covar[train_set,], tr = tx[train_set], yobs = Y_choice[train_set])

# estimate CATE for test set using train data
cate_esti_rf <- EstimateCate(xl_rf, covar[test_set,])
cate_esti_bart <- EstimateCate(xl_bart, covar[test_set,])

# estimate CATE for train set using test data
xl_rf2 <- X_RF(feat = covar[test_set,], tr = tx[test_set], yobs = Y_choice[test_set])
xl_bart2 <- X_BART(feat = covar[test_set,], tr = tx[test_set], yobs = Y_choice[test_set])

cate_esti_rf2 <- EstimateCate(xl_rf2, covar[train_set,])
cate_esti_bart2 <- EstimateCate(xl_bart2, covar[train_set,])

# combine both results into one vector
all_cate_esti_rf <- list()
all_cate_esti_rf[test_set] <- cate_esti_rf
all_cate_esti_rf[train_set] <- cate_esti_rf2
all_cate_esti_rf <- unlist(all_cate_esti_rf)

all_cate_esti_bart <- list()
all_cate_esti_bart[test_set] <- cate_esti_bart
all_cate_esti_bart[train_set] <- cate_esti_bart2
all_cate_esti_bart <- unlist(all_cate_esti_bart)

############################ GRF ########################################

# build causal forest using OOB samples
cf_oob <- causal_forest(
                    X = covar[train_set,],
                    Y = Y_choice[train_set],
                    W = tx[train_set],
                    num.trees = 6000

)
tau_hat_oob1 <- predict(cf_oob, estimate.variance = TRUE)
hist(tau_hat_oob1$predictions)

tau_hat_oob2 <- predict(cf_oob, covar[test_set,], estimate.variance = TRUE)
X_visualize <- "T_stage"
X_vis_id <- which(colnames(covar) %in% X_visualize)
plot(covar[test_set, X_vis_id], tau_hat_oob2$predictions, ylim = range(tau_hat_oob2$predictions, 0, 4.5), xlab = X_visualize, ylab = "tau", type = "l")
lines(covar[test_set, X_vis_id], pmax(0, covar[test_set, X_vis_id]), col = 2, lty = 2)

# estimation without splitting into test and train set
cf_noob <- causal_forest(
                    X = covar,
                    Y = Y_choice,
                    W = tx,
                    num.trees = 6000
)
tau_hat <- predict(cf_noob, estimate.variance = TRUE)

# it seems like the higher the T-stage, the more likely it's not going to have a positive treatment effect trend

ate_all <- average_treatment_effect(cf_noob, target.sample = "all")
ate_tx <- average_treatment_effect(cf_noob, target.sample = "treated")

tc <- test_calibration(cf_noob)

# Add confidence intervals for heterogeneous treatment effects; growing more trees is now recommended.
sigma_hat <- sqrt(tau_hat$variance.estimates)
CI_up <- mean(tau_hat$predictions) + sigma_hat * 1.96
CI_low <- mean(tau_hat$predictions) - sigma_hat * 1.96
summary(CI_up)
summary(CI_low)

############################# Powers methods ####################################

############ 1. PTO forest

p_forest1 <- PTOforest(x = covar[train_set,], tx = tx[train_set], y = Y_choice[train_set])
p_forest_pre1 <- predict(p_forest, newx = covar[test_set,])

p_forest2 <- PTOforest(x = covar[test_set,], tx = tx[test_set], y = Y_choice[test_set])
p_forest_pre2 <- predict(p_forest, newx = covar[train_set,])

all_p_forest_pre <- list()
all_p_forest_pre[test_set] <- p_forest_pre1
all_p_forest_pre[train_set] <- p_forest_pre2
all_p_forest_pre <- unlist(all_p_forest_pre)

############ 2. cross-validated causal boosting
cvcb1 <- cv.causalBoosting(x = covar[train_set,], tx = tx[train_set],, y = Y_choice[train_set],)
cvcb_pre1 <- predict(cvcb1, newx = covar[test_set,])

cvcb2 <- cv.causalBoosting(x = covar[train_set,], tx = tx[train_set],, y = Y_choice[train_set],)
cvcb_pre2 <- predict(cvcb2, newx = covar[test_set,])

all_cvcb_pre <- list()
all_cvcb_pre[test_set] <- cvcb_pre1
all_cvcb_pre[train_set] <- cvcb_pre2
all_cvcb_pre <- unlist(all_cvcb_pre)


############ 3. causal MARS

cm1 <- causalMARS(x = covar[train_set,], tx = tx[train_set],, y = Y_choice[train_set],)
cm_pre1 <- predict(cm1, newx = covar[test_set,])

cm2 <- causalMARS(x = covar[test_set,], tx = tx[test_set],, y = Y_choice[test_set],)
cm_pre2 <- predict(cm2, newx = covar[train_set,])

all_cm_pre <- list()
all_cm_pre[test_set] <- cm_pre1
all_cm_pre[train_set] <- cm_pre2
all_cm_pre <- unlist(all_cm_pre)

# compare all methods
vs_df <- data.frame()
vs_df <- rbind(vs_df, summary(all_cate_esti_rf))
vs_df <- rbind(vs_df, summary(all_cate_esti_bart))
vs_df <- rbind(vs_df, summary(tau_hat$predictions))
vs_df <- rbind(vs_df, summary(all_p_forest_pre))
vs_df <- rbind(vs_df, summary(all_cvcb_pre))
vs_df <- rbind(vs_df, summary(all_cm_pre))
colnames(vs_df) <- names(summary(all_cm_pre))

methods_name <- c("X_learner_RF", "X_learner_BART", "GRF", "PTO_forest", "CV_causal_boosting", "causal_MARS")
vs_df <- cbind(as.data.frame(methods_name), vs_df)