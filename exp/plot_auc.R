source("./bin/load_lib.R")
# SHAP values
suppressPackageStartupMessages({
library("SHAPforxgboost"); library("ggplot2"); library("xgboost")
library("data.table"); library("here"); library("svglite"); 
library(xgboostExplainer); library(rpart.plot) ; library(rattle); library(partykit)
})

min_mse_method <- read.csv("./res/best_tau_estimators.csv")
pos_report <- c("NCT00113763", "NCT00339183", "NCT00364013", "NCT00460265")
trial_choice <- filter(min_mse_method, trial %in% pos_report)
colnames(trial_choice)[1] <- "trialID"

# potentially pairing with NCT00364013
# for (j in 1:dim(trial_choice)[1]) {

    train_trial <- trial_choice[j,1]
    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))
    
    load(paste0("./bin/load_RCT/RCT_obj/", train_trial, ".RData"))
    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", train_trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", train_trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", train_trial, ".RData"))
    }
    X_train <- as.matrix(get(train_trial)[[1]])
    W_train <- get(train_trial)[[2]]

    outcome <- trial_choice[j,4]
    trial_best_method <- trial_choice[j,2]
    tau_train <- read.csv(paste0("./res/ite_tau_estimates/", train_trial, "_", outcome, "_", trial_best_method, "_tau_estimates.csv"))

    xgb_res <- try(readRDS(paste0("./dat/xgb_model/", train_trial, "_", outcome, "_", trial_best_method, "_xgb_model.rds"))) # generated from cvboost(x = X, y = tau[,1], objective="reg:squarederror")

    best_param <- xgb_res$best_param

    #
    # Load external validation data
    #
    test_trial <- "NCT00364013"
    load(paste0("./bin/load_RCT/RCT_obj/", test_trial, ".RData"))
    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", test_trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", test_trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", test_trial, ".RData"))
    }
    X_test <- as.matrix(get(test_trial)[[1]])
    W_test <- get(test_trial)[[2]]
    test_best_method <- "CBoost"
    tau_test <- read.csv(paste0("./res/ite_tau_estimates/", test_trial, "_", outcome, "_", test_best_method, "_tau_estimates.csv"))

    #
    # Make the two trials comparable
    #
    common_cov <- intersect(colnames(X_test), colnames(X_train))
    X_train <- X_train[,common_cov]
    X_test <- X_test[,common_cov]

    
    # Control
    train_ds <- data.frame(X_train, tau = as.factor(ifelse(tau_train[,1] > 0, "pos", "neg")))
    test_ds <- data.frame(X_test, tau = as.factor(ifelse(tau_test[,1] > 0, "pos", "neg")))
    rtree_model <- rpart(tau~., data=test_ds, control=rpart.control(maxdepth=5))

predictandCM<- function(amodel,data,modeltype)
{
  pred <-predict(amodel,data,type=modeltype)
  confusionMatrix(pred, reference=train_ds$tau,positive = "pos")
}
predictandCM(rtree_model,train_ds,"class")

predictandCM(rtree_model,test_ds,"class")
pred_test <- predict(rtree_model, test_ds ,"class")
confusionMatrix(pred_test, reference=test_ds$tau,positive = "pos")


train.control <- trainControl(
                           method = "repeatedcv",
                           number = 10, ## 10-fold CV
                           repeats = 3,## repeated three times
                           # USE AUC
                           summaryFunction = twoClassSummary, 
                           classProbs = TRUE
                           )

rpartFit1 <- train(tau~., data=train_ds,  
                 method = "rpart2", 
                 tuneLength = 10,
                 trControl = train.control,
                 metric = "ROC"
               )

rpartFit1

predictandCM(rpartFit1,test.data,"raw")

    pdf("test_plot.pdf")
    fancyRpartPlot(rtree_model) 
    dev.off()

    rpartTune <- train(tau~., data = train_ds, 
                    method = "rpart2", 
                    tuneLength = 100, 
        trControl = trainControl(method = "cv"))
    
    pdf("test_plot.pdf")
    plot(rpartTune)
    dev.off()

    rpartTree <- rpart(tau~., data=train_ds, maxdepth = 2)
    print(rpartTree)

    rpartTree2 <- as.party(rpartTree)
    pdf("test_plot.pdf")
    plot(rpartTree2) 
    dev.off()

    train.control <- trainControl(
                           method = "repeatedcv",
                           number = 10, ## 10-fold CV
                           repeats = 5,## repeated three times
                           # USE AUC
                           summaryFunction = defaultSummary, 
                           classProbs = TRUE
                           )

    reset.seed()
rpartFit1 <- train(tau~., data=train_ds, 
                 method = "rpart2", 
                 tuneLength = 6,
                 trControl = train.control,
                 metric = "RMSE"
               )
pdf("test_plot.pdf")
fancyRpartPlot(rpartFit1$finalModel)    
dev.off()

# } # trial loop