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

file_path <- "./dat/PDS/Colorec_Amgen_2005_262_NCT00115765/csv/"

datals <- c("corevar", "a_eendpt", "respeval","demo", "eligcrit", "lab", "lesion", "medhist", "chemotx", "radiotx", "surghist", "vitals_v")

for (sheet in datals) {
    assign(sheet, read.csv(paste0(file_path, sheet, ".csv")))
    if (sheet != "lab" & sheet != "lesion" & sheet != "medhist" & sheet != "vitals_v")
        assign(sheet, get(sheet) %>% replace_with_na_all(condition = ~.x %in% na_strings))
}

#
# Core variables all patients have data in
#
cv <- dplyr::select(corevar, all_of(c("SUBJID", "ACHEMM", "AGE", "SEXCD", "RACCAT", "B_LESNM", "B_LESNMT", "B_LESNMN", "B_METANM", "DURSURG", "PRADJYN", "B_LDHYN", "B_ECOG", "DIAGTYPE", "PDADEN", "PDPROMED", "KRAS")))

#
# Demographic variables
#
dem <- dplyr::select(demo, all_of(c("SUBJID", "B_WEIGHT", "B_HEIGHT", "B_BSA", "DIAGMONS", "METMONS", "CHILDPOT")))

#
# Baseline lab results
#
bl_lab <- filter(lab, LBBASE == "Y")
bl_lab <- dplyr::select(bl_lab, all_of(c("SUBJID", "LBTEST", "LBBLRES")))
bl_lab <- spread(bl_lab, LBTEST, LBBLRES)
bl_lab <- bl_lab %>% replace_with_na_all(condition = ~.x %in% na_strings)
colnames(bl_lab) <- gsub(" ", "_", colnames(bl_lab))

#
# Eligible patients
#
elig <- dplyr::select(eligcrit, all_of(c("SUBJID", "ELIGIBLE")))

# no need to include lesion data again

#
# Medical history (only number of unresolved body system entries included)
#
current_mh <- filter(medhist, MHSTAT == "Current or continuing")
current_mh <- current_mh %>% group_by(SUBJID) %>% add_count(MHTERM, name = "cond_count")
current_mh <- dplyr::select(current_mh, all_of(c("SUBJID", "MHTERM", "cond_count")))
current_mh <- unique(current_mh)
current_mh <- spread(current_mh, MHTERM, cond_count)
current_mh[is.na(current_mh)] <- 0
colnames(current_mh) <- gsub(" ", "_", colnames(current_mh))
colnames(current_mh)[2:dim(current_mh)[2]] <- paste0("MH_", colnames(current_mh)[2:dim(current_mh)[2]])

#
# Other prior cancer history, radiotherapy or chemo therapy histories
#
chany <- dplyr::select(chemotx, all_of(c("SUBJID", "CHANY5YR")))
chany <- unique(chany)

# radiotherapy history
raany <- dplyr::select(radiotx, all_of(c("SUBJID", "RAANY5YR")))
raany <- raany[-which(raany$SUBJID == "249196008" & raany$RAANY5YR == "N"),]
raany <- unique(raany)

# Surgery history
surgh <- surghist %>% group_by(SUBJID) %>% add_count(SXTYPE, name = "surgtype_count")
surgh <- dplyr::select(surgh, all_of(c("SUBJID", "SXTYPE", "surgtype_count")))
surgh <- unique(surgh)
surgh <- spread(surgh, SXTYPE, surgtype_count) %>% dplyr::select(-"<NA>")

#
# Vitals
#
vtl <- filter(vitals_v, VSBASE == "Y")
vtl <- vtl %>% replace_with_na_all(condition = ~.x %in% na_strings)
vtl <- vtl %>% group_by(SUBJID) %>% slice_max(DOSREFDY, n=1, with_ties = TRUE)
vtl <- dplyr::select(vtl, all_of(c("SUBJID", "VSTEST", "VSBLRES")))
vtl <- unique(vtl)
vtl <- spread(vtl, VSTEST, VSBLRES)
colnames(vtl) <- gsub(" ", "_", colnames(vtl))

#
# Join and combine all baseline characteristics
#
cleaned_dat <- c("bl_lab", "elig", "current_mh", "chany", "raany", "surgh", "vtl")

X <- left_join(cv, dem, by = c("SUBJID"))
for (datf in cleaned_dat) {
    X <- left_join(X, get(datf), by = c("SUBJID"))
    print(dim(unique(X)))
}
X <- unique(X)
sel_SUBJID <- X$SUBJID
X$SUBJID <- NULL
X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = as.data.frame(X), save_ddt = FALSE)

W <- as.numeric(corevar$ATRT == "panit. plus bevacizumab with chemotherapy")

#
# Export trial data
#
NCT00115765 <- list(X_imp, W)

# Objective response
endpt <- a_eendpt[a_eendpt$SUBJID %in% sel_SUBJID,]
OS <- data.frame(T = endpt$DTHMT, C = endpt$DTH) # primary outcome
PFS <- data.frame(T = endpt$PFSMTLR, C = endpt$PFSCR) # secondary outcome
TF <- data.frame(T = endpt$TFAILMT, C = endpt$TFAIL)

# Ordinal categorical outcome
ORR <- endpt$ROSLRCA # secondary outcome
rsps <- c("CR", "PR", "SD", "PD", "UE", "ND")
i <- 1
for (rsp_cat in rsps) {
    ORR[ORR == rsp_cat] <- i
    i <- i+1
}

NCT00115765_outcomes <- c("OS", "PFS", "TF", "ORR")

for (outcome in NCT00115765_outcomes) {
    if (outcome != "ORR") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = get(outcome)[,2], X = X_imp)))
    } else {
        ORR_Y_list <- list(ORR, NA, NA)

    }
}
