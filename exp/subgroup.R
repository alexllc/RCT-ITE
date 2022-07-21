source("./bin/load_lib.R")


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
        W <- get(trial)[[2]]
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

        res <- modav(resp = "y", trt = "treat", subgr = subgr, data = fitdat, covars = ~ x1 + x2, fitfunc = "lm")