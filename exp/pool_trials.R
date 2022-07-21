# Comparable trials pooled together
# mCRC with panitumumab: NCT00113763, NCT00115765, NCT00339183, NCT00364013, NCT00460265
# maybe (NCT00079274) but it uses Cetuximab rather than Panitumumab

adsl_trial <- c("NCT00113763", "NCT00364013")
long_trial <- c("NCT00339183", "NCT00115765")

for (trial in adsl_trial) {

    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))

    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }
    assign(paste0("X_", trial), as.matrix(get(trial)[[1]]))
    assign(paste0("W_", trial), get(trial)[[2]])
    assign(paste0("outcome_list_", trial), get(paste0(trial, "_outcomes")))

    for (outcome in get(paste0("outcome_list_", trial))) {
        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }

        # If at least one outcome in one trial is log transformed, all others should also be log transformed
        if(Y != "RSP") {
            Y <- log(Y)
            if (any(is.infinite(Y)))
                Y[which(is.infinite(Y))] <- 0
        }
        
        assign(paste0(trial, "_", outcome), Y)
    }
}

colnames(X_NCT00113763)[1] <- "PRSURG"

X_common <- intersect(colnames(X_NCT00113763), colnames(X_NCT00364013))
outcome_list_adsl <- intersect(outcome_list_NCT00113763, outcome_list_NCT00364013)

X <- rbind(dplyr::select(as.data.frame(X_NCT00113763), all_of(X_common)), dplyr::select(as.data.frame(X_NCT00364013), all_of(X_common)))
# LHD values are converted to categorical threshold
# Normal values of LHD are between 135 to 222U/L according to Mayo Clinic https://www.ncbi.nlm.nih.gov/books/NBK557536/
# Upper limit of normal range is 330 U/L for Cleveland Clinic https://pubmed.ncbi.nlm.nih.gov/17700600/
# We will use the Cleveland clinic upper limit here

# The ALP upper limit of normal range 147IU/L is also used here from Cleveland Clinic https://my.clevelandclinic.org/health/diagnostics/22029-alkaline-phosphatase-alp

# Creatinine For adult men, 0.74 to 1.35 mg/dL (65.4 to 119.3 micromoles/L); For adult women, 0.59 to 1.04 mg/dL (52.2 to 91.9 micromoles/L)

X <- X %>% mutate(LDH = Lactate_Dehydrogenase / 330, ALP = Alkaline_Phosphatase  / 147, Cre = ifelse(SEX == 2, Creatinine / 119.3, Creatinine / 91.9))

X <- dplyr::select(X, -c("Lactate_Dehydrogenase", "Alkaline_Phosphatase", "Creatinine"))

W <- c(W_NCT00113763, W_NCT00364013)

# Combine outcomes from two trials

for (outcome in get(paste0("outcome_list_", trial))) {
    assign(paste0(outcome, "_Y"), c())
    
    for (trial in adsl_trial) {
        Y <- get(paste0(trial, "_", outcome))

        assign(paste0(outcome, "_Y"), c(get(paste0(outcome, "_Y")), Y))

    }
}
    

for (outcome in outcome_list_adsl) {
    Y <- get(paste0(outcome, "_Y"))
    X <- as.matrix(X)
    rmv_id <- check_lm_rmv(X = X, Y = Y)

    X_rmvd <- X[-rmv_id,]
    Y_rmvd <- Y[-rmv_id]
    W_rmvd <- W[-rmv_id]

    tau <- perform_rstack(Y = Y_rmvd, X = X_rmvd, W = W_rmvd, trial_name = trial, outcome = outcome, imp_type_Y = imp_type_Y, tuned_cb_param = FALSE, tuned_cm_param = TRUE, tune_ptof_param = TRUE, perform_xb = FALSE, perform_cb = FALSE)
    
}


# For extracting lab test units

# for (trial in trial_adsl) {
#     source(paste0("./bin/load_RCT/load_", trial, ".R"))
#     lab_dt <- adlb %>% group_by(LBTEST) %>% slice(1)
#     lab_dt <- dplyr::select(lab_dt, all_of(c("LBTEST", "LBSTUNIT")))
#     write.csv(lab_dt, paste0("./dat/", trial, "_lab_dt.csv"), row.names = FALSE)
# }