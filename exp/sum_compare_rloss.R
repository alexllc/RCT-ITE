source("./bin/load_lib.R")

# Extract trial names from list of scripts
rloss_files <- list.files("./res/compare_rloss/")
trial_ls <- unique(unlist(lapply(strsplit(rloss_files, "_|\\."), function(X) X[1])))

rloss_list <- list()

for (trial in trial_ls) {
    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }

    trial_files <- rloss_files[grep(paste0("^", trial,"*"), rloss_files)]
    outcome_list <- unique(unlist(lapply(strsplit(trial_files, "_|\\."), function(X) X[5])))

    for (outcome in outcome_list) {
        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }

        iter_files <- trial_files[grep(paste0("_", outcome,"_"), trial_files)]
        iter_list <- unique(unlist(lapply(strsplit(iter_files, "_|\\."), function(X) X[6])))
        
        avg_rloss_list <- list()
        for (iter in iter_list) {
            rloss <- read.csv(paste0("./res/compare_rloss/", trial, "_model_sel_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"))
            avg_rloss_list <- append(avg_rloss_list, list(as.list(rloss$mse)))
        }
        rloss_iter_sum <- do.call(rbind.data.frame, avg_rloss_list)
        colnames(rloss_iter_sum) <- rloss$X

        rloss_list <- append(rloss_list, list(c(trial, outcome, as.list(colMeans(rloss_iter_sum, na.rm = TRUE)))))
        }
}

rloss_sum <- do.call(rbind.data.frame, rloss_list)
colnames(rloss_sum) <- c("trial", "outcome", rloss$X)

write.csv(rloss_sum, file = "./res/compare_rloss.csv", row.names = FALSE)