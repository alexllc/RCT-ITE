source("./bin/load_lib.R")
library(boot)

Q <- 4
min_mse_method <- read.csv("./res/best_tau_estimators.csv")
pos_report <- c("NCT00113763", "NCT00339183", "NCT00364013", "NCT00460265")
trial_choice <- filter(min_mse_method, mse < 150 & trial %in% pos_report)
hte_cstat <- read.csv("./cstat.csv")
colnames(hte_cstat) <- c("trial", "outcome", "c_benefit")
hte_cstat <- filter(hte_cstat, trial %in% trial_choice$trial)
colnames(trial_choice)[1] <- "trialID"

cstat_res <- data.frame(trial = character(), outcome = character(), cstat = numeric())

for (j in 9:dim(trial_choice)[1]) {
    
    trial <- trial_choice[j,1]
    message(paste0(rep("=", 80)))
    message(paste0("Loading: ", trial))
    message(paste0(rep("=", 80)))
    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    message("Data loaded.")
    
    outcome <- trial_choice[j,4]
    tau_method <- trial_choice[j,2]

    message(paste0("Outcome: ", outcome))

    # load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]
    Y_list <- get(paste0(outcome, "_Y_list"))

    imp_type_Y <- "efronYn"
    Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

    if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
        Y <- as.numeric(Y_list[[1]])
        imp_type_Y <- "efron"
    }

    time_indep <- X[,colnames(X)!="offtrt_reason"]
    cox_df <- cbind(time_indep, W)

    if (outcome == "RSP" | outcome == "ORR")
        next

    res.cox <- coxph(Surv(get(outcome)[["T"]], ceiling(get(outcome)[["C"]])) ~ . + W*(.), data = as.data.frame(cox_df))
    write.csv(tidy(res.cox), file = paste0("./res/cox_models/", trial, "_", outcome, "_", tau_method, "_cox_summary.csv"), row.names = FALSE)
}

    model <- glm(Y ~ .^2, data = as.data.frame(cox_df))

    con_pred <- predict(model)

    cstat <- calcualte_c_benefit(Y = Y, tau = con_pred, W = W)
    cstat_res <- rbind(cstat_res, c(trial, outcome, cstat))

    message("C for benefit is: ")
    print(cstat)
    print(cstat_res)
}


# Inverse probability censorship weighting


cens_prob <- ipwpoint(W, family = "binomial", link = "logit", numerator = NULL, denominator = formula(paste0("~", paste(colnames(cox_df)[1:10], collapse = "+"))), data = cox_df, trunc = NULL)


model <- glm(ceiling(time_event[["C"]]) ~ .^2, data = as.data.frame(cox_df), weights = cens_prob)

# Using Y prediction as CoxPH outcome
# [1] "NCT00052910"
# > outcome
# [1] "OS"
# > res
# [1] 0.7859413

# Using censoring status as CoxPH outcome
# > outcome
# [1] "OS"
# > res
0.2622074

#  0.472512
# trial
#  "NCT00364013"
#  outcome
#  "RSP"




    available_sample <- seq(1:dim(X)[1])
    holdout_sample <- sample(available_sample, floor(dim(X)[1] / Q))
    train_sample <- available_sample[!available_sample %in% holdout_sample]

    ho_X <- X[holdout_sample,]
    ho_W <- W[holdout_sample]
    ho_Y <- Y[holdout_sample]

    train_X <- X[train_sample,]
    train_W <- W[train_sample]
    train_Y <- Y[train_sample]
    
    
    results <- boot(data=mtcars, statistic=rsq, R=1000, formula=mpg~wt+disp)
