# Adapted for Tweedie likelihood from this very good post at https://www.simoncoulombe.com/2019/01/bayesian/
# objective function: we want to minimize the neg log-likelihood by tuning hyperparameters
obj.fun <- makeSingleObjectiveFunction(
  name = "causalboost_bayes",
  fn =   function(args){
    set.seed(42)
    cb <- causalBoosting(x = args["X"], 
                        tx = args["W"], 
                        y = args["Y"], 
                        num.trees = args["num.trees"], 
                        maxleaves = args["maxleaves"], 
                        eps = args["eps"],
                        splitSpread = args["splitSpread"])
    min_y_err <- cb$err.y[which(cb$err.y == min(cb$err.y))]
    tau_pred <- cb$tauhat[,500]
    
    rloss <- mean((args["Y"] - args["Y.hat"] - tau_pred * (args["W"] - args["W.hat"]))^2)

    return(rloss)
  },
  par.set = makeParamSet(
    makeNumericParam("num.trees",          lower = 500, upper = 100000),
    makeNumericParam("maxleaves",          lower = 1,     upper = 10),
    makeIntegerParam("eps",                lower= 0.01,      upper = 0.5),
    makeIntegerParam("splitSpread",       lower= 0.01,    upper = 0.5),
  ),
  minimize = TRUE ## rloss
)

#
# Optimization
#

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

