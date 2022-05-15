#' A sub-script of the main script load_NCT00115765.R to only select subjects receiving Oxaliplatin
#' Load PACCE trial results
#' Primary objective: 
#'  - Progression-Free Survival (Oxaliplatin)
#'  - Objective Tumor Response Through Week 12 (Irinotecan)
#' Secondary objective:
#'  - Overall Survival (Oxaliplatin)
#'  - Objective Tumor Response Rate (Oxaliplatin)
#'  - Time to Progression (Oxaliplatin)
#'  - Time to Treatment Failure (Oxaliplatin)
#'  - Overall Survival (Irinotecan)
#'  - Progression-free Survival (Irinotecan)
#'  - Objective Tumor Response Rate (Irinotecan)
#'  - Time to Progression (Irinotecan)
#'  - Time to Treatment Failure (Irinotecan)
#' 
#' 

source("./bin/load_RCT/load_NCT00115765.R")
oxa_subj <- which(corevar$CHEMO == "Oxaliplatin")

#
# Export trial data
#
NCT00115765_oxa <- list(NCT00115765[[1]][oxa_subj,], NCT00115765[[2]]W[oxa_subj])

# Objective response
endpt <- a_eendpt[oxa_subj,]
OS <- data.frame(T = endpt$DTHMT, C = endpt$DTH) # primary outcome
PFS <- data.frame(T = endpt$PFSMTCR, C = endpt$PFSCR) # secondary outcome
TF <- data.frame(T = endpt$TFAILMT, C = endpt$TFAIL)

# Ordinal categorical outcome
ORR <- endpt$ROSLRCA # secondary outcome
rsps <- c("CR", "PR", "SD", "PD", "UE", "ND")
i <- 1
for (rsp_cat in rsps) {
    ORR[ORR == rsp_cat] <- i
    i <- i+1
}

NCT00115765_oxa_outcomes <- c("OS", "PFS", "TF", "ORR")

for (outcome in NCT00115765_oxa_outcomes) {
    if (outcome != "ORR") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = get(outcome)[,2], X = X_imp)))
    } else {
        ORR_Y_list <- list(ORR, NA, NA)

    }
}
