# This is the software code used for generating the results in the final paper:
# Hohberg, M., Pütz, P., Kneib, T. (2020). Treatment effects beyond the mean 
# using distributional regression: Methods and guidance. Plos one 15.2 (2020):
# e0226514.

# clear workspace
remove(list = ls())

# choose seed, will be used to set seed below
my_seed <- 87

# set working directory where the data set is stored, the data set ("AER_credit.dta")
# can be obtained from https://www.aeaweb.org/articles?id=10.1257/aer.99.1.486, 
setwd("...")

# number of bootstrap replicates for pairs cluster bootstrap
boot_samples <- 499

# choose group of interest, default is the group of ineligibles,
# but you can also choose the eligibles treated or all people in
# treatment villages considered as treated
eligibles <- FALSE
all_considered <- FALSE

# function to load all required packages and install them if necessary
Install_And_Load <- function(Required_Packages) {
  Remaining_Packages <-
    Required_Packages[!(Required_Packages %in% installed.packages()[, "Package"])]
  
  
  if (length(Remaining_Packages)) {
    install.packages(Remaining_Packages)
  }
  for (package_name in Required_Packages)
  {
    library(package_name,
            character.only = TRUE,
            quietly = TRUE
    )
  }
}

# Specify the list of required packages to be installed and load
Required_Packages <-
  c(
    "foreign",
    "gamlss",
    "normalp",
    "ineq",
    "acid",
    "doParallel",
    "sandwich",
    "lmtest",
    "xtable",
    "doRNG"
  )

# Call the Function
Install_And_Load(Required_Packages)

# load data set
data <- read.dta("AER_credit.dta")
head(data)

# median income
median_inc <- median(data$Cind, na.rm = TRUE)

# poverty line: 60% of median
poverty_line <- median_inc * 0.6

# pairs cluster bootstrap function
clusbootreg <-
  ## input parameters
  # model: estimated model
  # -.formula: formula for the respective parameter of the response distribution
  # newdata_-: values of the explanatory variabels for control and treatment group to compute marginal treatment effect
  # start.from: vector of starting value for regression coefficients
  # family: response distribution
  # -.fix: is any of the four parameters assumed to fixed?
  # data: data set under consideration
  # cluster: cluster variable
  function(model,
           mu.formula = NULL,
           sigma.formula = NULL,
           nu.formula = NULL,
           tau.formula = NULL,
           newdata_control = NULL,
           newdata_treat = NULL,
           start.from = NULL,
           family = family,
           mu.fix = NULL,
           sigma.fix = NULL,
           nu.fix = NULL,
           tau.fix = NULL,
           nu.start = NULL,
           data,
           cluster) {
    # save cluster names
    clusters <- names(table(factor(cluster)))
    
    # set output to NA
    sterrs <- NA
    
    # do bootstrap until output is not NA anymore (i.e. draw boostrap samples until model with this sample converges)
    while (is.na(sterrs)[1]) {
      ## sample clusters
      # placeholder for cluster indicator
      boot.ind <- c()
      # sample index
      index <-
        sample(1:length(clusters), length(clusters), replace = TRUE)
      # loop through number of clusters
      for (k in 1:length(clusters)) {
        # save indicator of selected cluster
        obs.sel <- which(cluster == unique(cluster)[index[k]])
        # append the indicator to the other indicators
        boot.ind <- c(boot.ind, obs.sel)
      }
      # finally: bootstrap data set is build from sampled clustered (according to indicator) from original dataset
      bootdat <- data[boot.ind, ]
      
      # estimate gamlss
      err <- tryCatch({
        modgamlss <-
          (
            gamlss(
              formula = mu.formula,
              sigma.formula = sigma.formula,
              nu.formula = nu.formula,
              tau.formula = tau.formula,
              trace = FALSE,
              family = family,
              # estimation algorithm details
              method = mixed(10, 80),
              control = gamlss.control(
                mu.step = 0.1,
                sigma.step = 0.1,
                tau.step = 0.1,
                c.crit = 0.01
              ),
              data = bootdat
            )
          )
      }, error = function(e)
        e)
      # if convergence failure, save NA for this bootstrap iteration
      if (inherits(err, "error")) {
        sterrs <- NA
      } else {
        if (modgamlss$converged == TRUE) {
          ## marginal treatment effects at means, have to be adjusted to the
          ## respective response distribution and the respective
          ## distributional measures of interest, here: Burr distribution (GB2)
          ## pred_"group"_"parameter" indicates predicted respective parameter for t or c group
          pred_control_mu <-
            predict(
              modgamlss,
              newdata = newdata_control,
              what = "mu",
              type = "response",
              data = bootdat
            )
          pred_control_sigma <-
            predict(
              modgamlss,
              newdata = newdata_control,
              what = "sigma",
              type = "response",
              data = bootdat
            )
          pred_control_tau <-
            predict(
              modgamlss,
              newdata = newdata_control,
              what = "tau",
              type = "response",
              data = bootdat
            )
          pred_treat_mu <-
            predict(
              modgamlss,
              newdata = newdata_treat,
              what = "mu",
              type = "response",
              data = bootdat
            )
          pred_treat_sigma <-
            predict(
              modgamlss,
              newdata = newdata_treat,
              what = "sigma",
              type = "response",
              data = bootdat
            )
          pred_treat_tau <-
            predict(
              modgamlss,
              newdata = newdata_treat,
              what = "tau",
              type = "response",
              data = bootdat
            )
          # sequence of numbers for which should be predicted
          x <- seq(0.0001, 10000, by = 0.1)
          # compute densities of Burr distribution for control and treatment groups
          dens1 <-
            dGB2(
              x,
              mu = (pred_control_mu),
              sigma = pred_control_sigma,
              nu = 1,
              tau = pred_control_tau
            )
          dens2 <-
            dGB2(
              x,
              mu = (pred_treat_mu),
              sigma = pred_treat_sigma,
              nu = 1,
              tau = pred_treat_tau
            )
          # compute measures of interest (here: gini, atkinson etc.): marginal effects at means
          x <- seq(0.000001, 10000, by = 0.1)
          sterrs <-
            c(
              gini.den(incs = x, dens = dens1)$Gini,
              gini.den(incs = x, dens = dens2)$Gini,
              arithmean.GB2(
                pred_control_mu,
                pred_control_sigma,
                1,
                pred_control_tau
              ),
              arithmean.GB2(pred_treat_mu, pred_treat_sigma, 1, pred_treat_tau),
              sd.GB2(
                pred_control_mu,
                pred_control_sigma,
                1,
                pred_control_tau
              )^2,
              sd.GB2(pred_treat_mu, pred_treat_sigma, 1, pred_treat_tau)^
                2,
              atkinson.GB2(
                pred_control_mu,
                pred_control_sigma,
                1,
                pred_control_tau,
                epsilon = 1
              ),
              atkinson.GB2(
                pred_treat_mu,
                pred_treat_sigma,
                1,
                pred_treat_tau,
                epsilon = 1
              ),
              atkinson.GB2(
                pred_control_mu,
                pred_control_sigma,
                1,
                pred_control_tau,
                epsilon = 2
              ),
              atkinson.GB2(
                pred_treat_mu,
                pred_treat_sigma,
                1,
                pred_treat_tau,
                epsilon = 2
              ),
              entropy.GB2(
                pred_control_mu,
                pred_control_sigma,
                1,
                pred_control_tau,
                alpha = 1
              ),
              entropy.GB2(
                pred_treat_mu,
                pred_treat_sigma,
                1,
                pred_treat_tau,
                alpha = 1
              ),
              pGB2(
                poverty_line,
                pred_control_mu,
                pred_control_sigma,
                1,
                pred_control_tau
              ),
              pGB2(
                poverty_line,
                pred_treat_mu,
                pred_treat_sigma,
                1,
                pred_treat_tau
              )
            )
        } else {
          sterrs <- NA
        }
      }
    }
    sterrs
  }

###################
# determine time point of survey: Restrict data set to November 1999
dat <- subset(data, t == 10)

# remove consumption values not smaller than 10000
da <- subset(dat, Cind < 10000)

# if eligibles considered (see above)
if (eligibles == TRUE) {
  da$treatnp <- da$treatp
}

# if all people in treatment villages are considered treated (see above)
if (all_considered == TRUE) {
  da$treat <- rep(NA, dim(da)[1])
  da$treat[da$treatnp == 1] <- 1
  da$treat[da$treatp == 1] <- 1
  da$treat[da$treatnp == 0] <- 0
  da$treat[da$treatp == 0] <- 0
  da$treatnp <- da$treat
}

# retain variables of interest
da <-
  da[, c(
    "Cind",
    "treatnp",
    "indice",
    "hectareas",
    "hhhsex",
    "hhhage",
    "hhhalpha",
    "vhhnum",
    "yycali_1",
    "p16",
    "village"
  )]

# remove missing values
da <- na.omit(da)

# clear workspace
rm(dat,data)

# center variables
da$yycali_1 <- da$yycali_1 - mean(da$yycali_1)
da$vhhnum <- da$vhhnum - mean(da$vhhnum)
da$hhhage <- da$hhhage - mean(da$hhhage)
da$hhhsex <- as.factor(da$hhhsex)

# delete seldom/uninformative answers
da$hhhsex[da$hhhsex == 9] <- NA
da$p16[da$p16 == "nr"] <- NA
da$hhhalpha[da$hhhalpha == "nr"] <- NA

# delete NAs
da <- na.omit(da)

# recode pecularities
da$hhhalpha <- ifelse(da$hhhalpha == "no", "no", "si")
da$p16 <- ifelse(da$p16 == "no", "no", "si")

# how many consumption values equal to zero?
sum(da$Cind == 0)

# remove few observations with zero consumption for making log transformation possible,
# there might be also subject-matter reasons to remove these values
da <- subset(da, Cind > 0)

# estimation sample ginis in treatment and control
gini(da$Cind[da$treatnp == 0])$Gini
gini(da$Cind[da$treatnp == 1])$Gini

## replication of original model by Angelucci
form <-
  as.formula(Cind ~ treatnp + indice + hectareas + hhhsex + hhhage + vhhnum +
               yycali_1 + p16)

# coefficients
pm1 <- lm(form, data = da)
summary(pm1)

# cluster-robust inference:
G <- length(unique(da$village))
N <- length(da$village)
dfa <- (G / (G - 1)) * (N - 1) / pm1$df.residual
village_c_vcov <-
  dfa * vcovHC(pm1,
               type = "HC0",
               cluster = "village",
               adjust = T
  )
coeftest(pm1, vcov = village_c_vcov)

### FIGURE 1 ###
# describe dependent variable
summary(log(da$Cind))
histDist(log(da$Cind), "NO")
par(mfrow = c(1, 2))
hist((da$Cind),
     breaks = 30,
     freq = FALSE,
     xlim = c(0, 2000),
     xlab = "Consumption",
     main = "",
     col = "#d4d4f7",
     yaxt = "n"
)
axis(2,
     at = c(0, 0.0010, 0.0020, 0.0030),
     labels = c(0, 0.001, 0.002, 0.003)
)
hist(
  log(da$Cind),
  breaks = 30,
  xlim = c(2, 8),
  freq = FALSE,
  xlab = "Log Consumption",
  main = "",
  col = "#d4d4f7"
)
x <- seq(0, 10, length = 10000)
mu <- mean(log(da$Cind))
sigma <- sd(log(da$Cind))
points(x,
       dnorm(x, mu, sigma),
       type = "l",
       lwd = 2,
       col = "#df4620"
)

# further investigation of dependent variable grouped by treatment variable
par(mfrow=c(1,1))
plot(density(log(da$Cind[da$treatnp == 0])))
lines(density(log(da$Cind[da$treatnp == 1])), col = 2)

# save village variable
village <- da$village
# delete village variable from data frame
da <- da[, !(names(da) %in% "village")]

# log-normal model
modgama <- gamlss(Cind ~ .,
                  sigma.fo = ~.,
                  family = LOGNO,
                  data = da
)

# model diagnostics
plot(
  modgama,
  parameters = par(
    mfrow = c(2, 2),
    mar = par("mar") + c(0, 1, 0, 0),
    col.axis = "black",
    col.main = "black",
    col.lab = "black",
    col = "black",
    bg = "white",
    bty = "l",
    cex.axis = 1
  )
)

# burr distribution
modgam2 <-
  gamlss(
    Cind ~ .,
    sigma.fo = ~.,
    ~1,
    ~.,
    family = GB2(mu.link = "log"),
    nu.fix = TRUE,
    nu.start = 1,
    method = mixed(2, 20),
    start.from = modgama,
    data = da
  )
# model diagnostics
plot(
  modgam2,
  parameters = par(
    mfrow = c(2, 2),
    mar = par("mar") + c(0, 1, 0, 0),
    col.axis = "black",
    col.main = "black",
    col.lab = "black",
    col = "black",
    bg = "white",
    bty = "l",
    cex.axis = 1
  )
)


### FIGURE 2 (only the plot in the respective lower right panel) ###
# store residuals
residlogno <- resid(modgama)
residburr <- resid(modgam2)
# plot qq-plots of both models
par(mfrow=c(1,2))
qqnorm(residlogno, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
       col = "darkgreen")
lines(residlogno, residlogno, col = "red", lwd = 0.4, cex = 0.4)
qqnorm(residburr, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles", plot.it = TRUE, frame.plot = TRUE, 
       col = "darkgreen")
lines(residburr, residburr, col = "red", lwd = 0.4, cex = 0.4)


### TABLE 1 ###
# summary of quantile residuals for both models
table_logno_res <- plot(modgama)
table_burr_res <- plot(modgam2)

## predict the conditional distributions at mean values for treatment and control group
# determine values of variable for control group
newdata_control <-
  data.frame(
    treatnp = 0,
    indice = mean(da$indice),
    hectareas = mean(da$hectareas),
    hhhsex = 1,
    hhhage = mean(da$hhhage),
    hhhalpha = "si",
    vhhnum = mean(da$vhhnum),
    yycali_1 = mean(da$yycali_1),
    p16 = "no"
  )
# determine values of variable for treatment group
newdata_treat <-
  data.frame(
    treatnp = 1,
    indice = mean(da$indice),
    hectareas = mean(da$hectareas),
    hhhsex = 1,
    hhhage = mean(da$hhhage),
    hhhalpha = "si",
    vhhnum = mean(da$vhhnum),
    yycali_1 = mean(da$yycali_1),
    p16 = "no"
  )

## predict single parameters of burr distribution with these values
# for control group
pred_control_mu <-
  predict(modgam2,
          newdata = newdata_control,
          what = "mu",
          type = "response"
  )
pred_control_sigma <-
  predict(modgam2,
          newdata = newdata_control,
          what = "sigma",
          type = "response"
  )
pred_control_tau <-
  predict(modgam2,
          newdata = newdata_control,
          what = "tau",
          type = "response"
  )
# for treatment group
pred_treat_mu <-
  predict(modgam2,
          newdata = newdata_treat,
          what = "mu",
          type = "response"
  )
pred_treat_sigma <-
  predict(modgam2,
          newdata = newdata_treat,
          what = "sigma",
          type = "response"
  )
pred_treat_tau <-
  predict(modgam2,
          newdata = newdata_treat,
          what = "tau",
          type = "response"
  )

### FIGURE 3 ###
# plot estimated distributions
par(mfrow = c(1, 1))
x <- (seq(0, 600, length = 1000))
plot(
  x,
  dGB2(
    x,
    mu = (pred_control_mu),
    sigma = pred_control_sigma,
    nu = 1,
    tau = pred_control_tau
  ),
  type = "l",
  ylab = "Density",
  xlab = "Food consumption",
  col = "black",
  #ylim = c(0, 0.007),
  bty = "l",
  yaxt = "n"
)
axis(2,
     at = c(0, 0.002, 0.004, 0.006, 0.008),
     labels = c(0, 0.002, 0.004, 0.006, 0.008))
points(x,
       dGB2(
         x,
         mu = (pred_treat_mu),
         sigma = pred_treat_sigma,
         nu = 1,
         tau = pred_treat_tau
       ),
       col= "red",
       type = "l")
abline(v=poverty_line, col = "blue", lty = 2)
legend(
  330,
  0.005,
  c("control villages", "treatment villages", "poverty line"),
  col = c(1, "red", "blue"),
  lty = c(1,1,2),
  border = NULL
)

## compute point estimates of interest for control and treatment group
# sequence of x values
x <- seq(0.000001, 10000, by = 0.1)
# densities for these x values
dens1 <-
  dGB2(
    x,
    mu = (pred_control_mu),
    sigma = pred_control_sigma,
    nu = 1,
    tau = pred_control_tau
  )
dens2 <-
  dGB2(
    x,
    mu = (pred_treat_mu),
    sigma = pred_treat_sigma,
    nu = 1,
    tau = pred_treat_tau
  )
# compute point estimates for inequality measures
point_estimates_per_group <-
  c(
    gini.den(incs = x, dens = dens1)$Gini,
    gini.den(incs = x, dens = dens2)$Gini,
    arithmean.GB2(
      pred_control_mu,
      pred_control_sigma,
      1,
      pred_control_tau
    ),
    arithmean.GB2(pred_treat_mu, pred_treat_sigma, 1, pred_treat_tau),
    sd.GB2(
      pred_control_mu,
      pred_control_sigma,
      1,
      pred_control_tau
    )^2,
    sd.GB2(pred_treat_mu, pred_treat_sigma, 1, pred_treat_tau)^
      2,
    atkinson.GB2(
      pred_control_mu,
      pred_control_sigma,
      1,
      pred_control_tau,
      epsilon = 1
    ),
    atkinson.GB2(
      pred_treat_mu,
      pred_treat_sigma,
      1,
      pred_treat_tau,
      epsilon = 1
    ),
    atkinson.GB2(
      pred_control_mu,
      pred_control_sigma,
      1,
      pred_control_tau,
      epsilon = 2
    ),
    atkinson.GB2(
      pred_treat_mu,
      pred_treat_sigma,
      1,
      pred_treat_tau,
      epsilon = 2
    ),
    entropy.GB2(
      pred_control_mu,
      pred_control_sigma,
      1,
      pred_control_tau,
      alpha = 1
    ),
    entropy.GB2(
      pred_treat_mu,
      pred_treat_sigma,
      1,
      pred_treat_tau,
      alpha = 1
    ),
    pGB2(
      poverty_line,
      pred_control_mu,
      pred_control_sigma,
      1,
      pred_control_tau
    ),
    pGB2(
      poverty_line,
      pred_treat_mu,
      pred_treat_sigma,
      1,
      pred_treat_tau
    )
  )

# remove some objects from workspace
rm(modgama)

# parallel computating for obtaining pairs bootstrap standard errors for marginal treatment effects at means
# how many cores to use
n.cores <- detectCores() - 1
# initializing the parallel system
cl <- makeCluster(n.cores)
registerDoParallel(cl)

# deliver all relevant parameters/arguments (see above)
output.parallel <-
  foreach(
    z = 1:boot_samples,
    .combine = "rbind",
    .packages = c("gamlss", "acid")
  ) %dorng% {
    set.seed(my_seed+z)
    print(z)
    variances <-
      clusbootreg(
        model = modgam2,
        mu.formula = Cind ~ .,
        sigma.formula = ~.,
        nu.formula = ~1,
        tau.formula = ~.,
        data = da,
        cluster = village,
        newdata_control = newdata_control,
        newdata_treat = newdata_treat,
        nu.fix = TRUE,
        family = GB2
      )
  }
stopCluster(cl)
# remove NAs from output (non-converging iterations)
output.parallel <- na.omit(output.parallel)
# note: not all models have necessarily converged, how many have converged?
print(c("number of valid bootstrap samples:", dim(output.parallel)[1]))

# check bootstrap statistics for symmetry / outliers
par(mfrow = c(2, (dim(output.parallel)[2]) / 2))
apply(output.parallel, 2, summary)
effect <-
  c(
    "Gini",
    "Mean",
    "Variance",
    "Atkinson (e=1)",
    "Atkinson (e=2)",
    "Theil",
    "Vulnerability"
  )

par(mfrow = c(2, 1), bg = "white")
for (j in 1:((dim(output.parallel)[2]) / 2))
{
  boxplot(output.parallel[, ((j - 1) * 2 + 2)] - output.parallel[, ((j - 1) *
                                                                      2 + 1)],
          horizontal = TRUE,
          main = effect[j]
  )
  hist(
    output.parallel[, ((j - 1) * 2 + 2)] - output.parallel[, ((j - 1) * 2 +
                                                                1)],
    breaks = 20,
    xlab = "",
    main = effect[j],
    col = "#d4d4f7"
  )
}
### FIGURE A can be produced by the following command if the ineligbles are consideres (see above)  ###
par(mfrow = c(2, 1), bg = "white", mai = c(0.5, 1.0, 0.2, 0.5))
boxplot(output.parallel[, 2] - output.parallel[, 1],
        horizontal = TRUE,
        #   xaxt = "n",
        main = ""
)
# axis(1, c(-0.02,0.00,0.02,0.04,0.05))
hist(
  output.parallel[, 2] - output.parallel[, 1],
  breaks = "Scott",
  main = "",
  xlab = "",
  xaxt = "n",
  col = "#d4d4f7"
)
axis(1, c(-0.02, 0.00, 0.02, 0.04))

# check percentile interval boundaries for convergence
# remove first burn_in observations to obtain less noisy trajectories
# works only if burn_in is larger than number of valid bootstrap estimates
burn_in <- 10
par(mfrow = c(2, dim(output.parallel)[2] / 2))
for (j in 1:((dim(output.parallel)[2]) / 2))
{
  rol_lb <- numeric(dim(output.parallel)[1])
  rol_ub <- numeric(dim(output.parallel)[1])
  for (i in 1:dim(output.parallel)[1])
  {
    rol_lb[i] <-
      quantile(output.parallel[1:i, ((j - 1) * 2 + 1)] - output.parallel[1:i, ((j -
                                                                                  1) * 2 + 2)], probs = 0.025)
    rol_ub[i] <-
      quantile(output.parallel[1:i, ((j - 1) * 2 + 1)] - output.parallel[1:i, ((j -
                                                                                  1) * 2 + 2)], probs = 0.975)
  }
  plot((burn_in + 1):dim(output.parallel)[1],
       rol_lb[-(1:burn_in)],
       type = "l",
       xlab = "",
       ylab = "",
       col = j,
       main = effect[j],
       bty = "l"
  )
  plot((burn_in + 1):dim(output.parallel)[1],
       rol_ub[-(1:burn_in)],
       type = "l",
       xlab = "",
       ylab = "",
       col = j,
       main = effect[j],
       bty = "l"
  )
}
### FIGURE B can be produced by the following command if the ineligbles are consideres (see above)  ###
par(mfrow = c(1, 2), bg = "white", mai = c(1.0, 1.0, 0.5, 0.5))

rol_lb <- numeric(dim(output.parallel)[1])
rol_ub <- numeric(dim(output.parallel)[1])
for (i in 1:dim(output.parallel)[1])
{
  rol_lb[i] <-
    quantile(output.parallel[1:i, 1] - output.parallel[1:i, 2], probs = 0.025)
  rol_ub[i] <-
    quantile(output.parallel[1:i, 1] - output.parallel[1:i, 2 ], probs = 0.975)
}
plot((burn_in + 1):dim(output.parallel)[1],
     rol_lb[-(1:burn_in)],
     type = "l",
     xlab = "Bootstrap replicates",
     ylab = "Mean estimate",
     col = 1,
     main = "Lower bound",
     bty = "l"
)
plot((burn_in + 1):dim(output.parallel)[1],
     rol_ub[-(1:burn_in)],
     type = "l",
     xlab = "Bootstrap replicates",
     ylab = "Mean estimate",
     col = 1,
     main = "Upper bound",
     bty = "l"
)



# marginal effects at means on different measures: difference of point estimates for treatment and control group:
le <- length(point_estimates_per_group)
measures <- numeric(le / 2)
for (i in 1:(le / 2))
{
  measures[i] <- point_estimates_per_group[i * 2] - point_estimates_per_group[i * 2 - 1]
}

### TABLE 2, 3 or 4; dependent on which group is the treatment group:
### eligibles==FALSE and all_cosidered==FALSE gives TABLE 2,
### eligibles==TRUE and all_cosidered==FALSE gives TABLE 3 and
### all_cosidered=TRUE gives TABLE 4
## print point estimates and 95% percentile confidence intervals
conf_inel <- data.frame(matrix(0, length(measures), 3, byrow = TRUE))
for (j in 1:((dim(output.parallel)[2]) / 2))
{
  conf_inel[j, ] <-
    c(
      measures[j],
      quantile(output.parallel[, ((j - 1) * 2 + 2)] - output.parallel[, ((j -
                                                                            1) * 2 + 1)], probs = c(0.025, 0.975))
    )
}
rownames(conf_inel) <- c(
  "MTE on Gini coefficient",
  "MTE on mean",
  "MTE on variance",
  "MTE on Atkinson index (e=1)",
  "MTE on Atkinson index (e=2)",
  "MTE on Theil index",
  "MTE on vulnerability"
)
colnames(conf_inel) <- c("Estimate", "Lower Bound", "Upper Bound")

latextable <-
  xtable(print(round(conf_inel[c(2, 3, 1, 4:7), ], digits = 3), row.names = TRUE),
         digits = 3
  )
align(latextable) <- xalign(latextable)
latextable
