source("./bin/load_lib.R")


# Remove covariates with high VIF values
check_colinearity <- function(Y = NULL, X = NULL) { # using 5 as threshold at the moment https://www.investopedia.com/terms/v/variance-inflation-factor.asp
    lm_df <- data.frame(outcome = Y, X)
    model <- lm(outcome~., data = lm_df)
    predictions <- model %>% predict(lm_df)
    print(data.frame(RMSE = RMSE(predictions, lm_df$outcome), R2 = caret::R2(predictions, lm_df$outcome)))

    vif_res <- car::vif(model)
    return(vif_res)
}


# Extract trial names from list of scripts
trial_scripts <- list.files("./bin/load_RCT")
trial_ls <- trial_scripts[grep("^load*", trial_scripts)]
trial_ls <- unlist(lapply(strsplit(trial_ls, "_|\\."), function(X) X[2]))

for (trial in trial_ls) {
    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))

    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }
    X <- as.matrix(get(trial)[[1]])

    res <- cor(X)
    pdf(paste0("./res/corr_plots/", trial, "_corr_plot.pdf"))
    # Insignificant correlation are crossed
    corrplot(res, type="upper", order="hclust", 
            tl.cex = 0.5)

    dev.off()
    outcome_list <- get(paste0(trial, "_outcomes"))
    for (outcome in outcome_list) {
        message(paste0("Processing outcome: ", outcome))

        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }        
        
        vif_res <- check_colinearity(Y = Y, X = X)
        print(vif_res)
    }
}