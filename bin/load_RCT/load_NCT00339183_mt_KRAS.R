source("./bin/load_RCT/load_NCT00339183.R")

mtkras_sel <- corevar$KRAS == "Mutant"
mtkras_SUB <- corevar$SUBJID[mtkras_sel]
#
# Export trial data
#
NCT00339183_mt_KRAS <- list(X_imp[mtkras_sel,], W[mtkras_sel])

# Objective response
endpt <- a_eendpt[a_eendpt$SUBJID %in% mtkras_SUB,]
endpt <- left_join(endpt, disposit, by = c("SUBJID"))
PFS <- data.frame(T = endpt$PDDYLR / 30.4167, C = endpt$PDLR) # secondary outcome
OS <- data.frame(T = endpt$DTHDY / 30.4167 , C = endpt$DTH) # primary outcome

# one DTHDY entry was missing, last contact day relative to 1st dose will be used instead
OS[which(is.na(OS$T)), 1] <- endpt$LASTCTDY[which(is.na(OS$T))] / 30.4167 # convert days to months


# Ordinal categorical outcome
ORR <- endpt$ROSLRCA # secondary outcome
rsps <- c("CR", "PR", "SD", "PD", "UE", "ND")
i <- 1
for (rsp_cat in rsps) {
    ORR[ORR == rsp_cat] <- i
    i <- i+1
}

NCT00339183_mt_KRAS_outcomes <- c("OS", "PFS", "ORR")

for (outcome in NCT00339183_mt_KRAS_outcomes) {
    if (outcome != "ORR") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = get(outcome)[,2], X = X_imp)))
    } else {
        ORR_Y_list <- list(ORR, NA, NA)

    }
}