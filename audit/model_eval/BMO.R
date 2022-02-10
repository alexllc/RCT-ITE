library(CASdatasets)
library(dplyr)
library(tibble)
library(magrittr)
library(ggplot2)
library(scatterplot3d)
library(kableExtra)
library(tidyr)
library(mlrMBO)
library(ParamHelpers)
library(DiceKriging)
library(smoof)
library(xgboost)

data("freMPL1")
data("freMPL2")
data("freMPL3")
fre_df <- rbind(freMPL1, freMPL2, freMPL3 %>% select(-DeducType))
rm(freMPL1, freMPL2, freMPL3)

gridExtra::grid.arrange(
  fre_df %>%
    filter(ClaimAmount > 0) %>%
    ggplot(aes(x = ClaimAmount)) +
    geom_density() +
    ggtitle("Observed Loss Distribution"),

  fre_df %>%
    filter(ClaimAmount > 0, ClaimAmount < 1.5e4) %>%
    ggplot(aes(x = ClaimAmount)) +
    geom_density() +
    ggtitle("Observed Severity Distribution"),
  nrow = 1
)

min(fre_df$ClaimAmount)

## [1] -3407.7

sum(fre_df$ClaimAmount < 0)

## [1] 690

fre_df <- fre_df %>%
  mutate(ClaimAmount = case_when(ClaimAmount < 0 ~ 0, TRUE ~ ClaimAmount)) %>%
  mutate(VehMaxSpeed_num = sub(".*-", "", VehMaxSpeed) %>% substr(., 1, 3)%>% as.numeric,
         VehAge_num = sub("*.-", "", VehAge) %>% sub('\\+', '', .) %>% as.numeric,
         VehPrice_num = as.integer(VehPrice)) %>% # The factor levels appear to be ordered so I will use this
  group_by(SocioCateg) %>% # high cardinality, will encode as a proportion of total
  mutate(SocioCateg_prop =  (sum(n()) / 4) / nrow(.) * 1e5) %>% 
  ungroup()

## matrices, no intercept needed and don't forget to exclude post-dictors
fre_mat <- model.matrix(ClaimAmount ~ . -1 -ClaimInd -Exposure -RecordBeg 
                        -RecordEnd - VehMaxSpeed -VehPrice -VehAge -SocioCateg,
                        data = fre_df)
## xgb.DMatrix, faster sparse matrix
fre_dm <- xgb.DMatrix(data = fre_mat, 
                      label = sample(fre_df$ClaimAmount, dim(fre_mat)[1]), 
                      base_margin = log(fre_df$Exposure)) ## base-margin == offset
                                                          ## we use log earned exposure because the xgboost Tweedie
                                                          ## implementation includes a log-link for the variance power

#
# Building objective function
#

# Adapted for Tweedie likelihood from this very good post at https://www.simoncoulombe.com/2019/01/bayesian/
# objective function: we want to minimize the neg log-likelihood by tuning hyperparameters
obj.fun <- makeSingleObjectiveFunction(
  name = "xgb_cv_bayes",
  fn =   function(x){
    set.seed(42)
    cv <- xgb.cv(params = list(
      booster          = "gbtree",
      eta              = x["eta"],
      max_depth        = x["max_depth"],
      min_child_weight = x["min_child_weight"],
      gamma            = x["gamma"],
      subsample        = x["subsample"],
      colsample_bytree = x["colsample_bytree"],
      max_delta_step   = x["max_delta_step"],
      tweedie_variance_power = x["tweedie_variance_power"],
      objective        = 'reg:tweedie', 
      eval_metric     = paste0("tweedie-nloglik@", x["tweedie_variance_power"])),
      data = dm, ## must set in global.Env()
      nround = 200, ## Set this large and use early stopping
      nthread = 26, ## Adjust based on your machine
      nfold =  5,
      prediction = FALSE,
      showsd = TRUE,
      early_stopping_rounds = 25, ## If evaluation metric does not improve on out-of-fold sample for 25 rounds, stop
      verbose = 1,
      print_every_n = 500)

    cv$evaluation_log %>% pull(4) %>% min  ## column 4 is the eval metric here, tweedie negative log-likelihood
  },
  par.set = makeParamSet(
    makeNumericParam("eta",                    lower = 0.005, upper = 0.01),
    makeNumericParam("gamma",                  lower = 1,     upper = 5),
    makeIntegerParam("max_depth",              lower= 2,      upper = 10),
    makeIntegerParam("min_child_weight",       lower= 300,    upper = 2000),
    makeNumericParam("subsample",              lower = 0.20,  upper = .8),
    makeNumericParam("colsample_bytree",       lower = 0.20,  upper = .8),
    makeNumericParam("max_delta_step",         lower = 0,     upper = 5),
    makeNumericParam("tweedie_variance_power", lower = 1.75,   upper = 1.85)
  ),
  minimize = TRUE ## negative log likelihood
)

do_bayes <- function(n_design = NULL, opt_steps = NULL, of = obj.fun, seed = 42) {
  set.seed(seed)

  des <- generateDesign(n=n_design,
                        par.set = getParamSet(of),
                        fun = lhs::randomLHS)

  control <- makeMBOControl() %>%
    setMBOControlTermination(., iters = opt_steps)

  ## kriging with a matern(3,2) covariance function is the default surrogate model for numerical domains
  ## but if you wanted to override this you could modify the makeLearner() call below to define your own
  ## GP surrogate model with more or lesss smoothness, or use an entirely different method
  run <- mbo(fun = of,
             design = des,
             learner = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2", control = list(trace = FALSE)),
             control = control, 
             show.info = TRUE)

}

dm <- fre_dm
runs <- do_bayes(n_design = 15, of = obj.fun, opt_steps = 10, seed = 42)