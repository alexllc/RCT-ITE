#' Script to load the SPECTRUM trial NCT00460265
#' Primary obj: overall survival (DTHDY, DTH)
#' Secondary obj: 
#' - Overall Objective Response Rate, (RDTLR, )
#' - Duration of Response, Time to Progression,
#' - Time to Response, 
#' - Progression Free Survival


file_path <- "./dat/PDS/HeadNe_Amgen_2007_265_NCT00460265/csv/"

datals <- c("corevar", "aeendpt", "demo", "lab", "eligcrit", "lesion", "chemotx", "othcantx", "surghist", "ae")

# set to NA in R language for any type of missing or not otherwise specified entries
# na strings loaded in the load_lib.R script
for (sheet in datals) {
    assign(sheet, read.csv(paste0(file_path, sheet, ".csv")))
    if (sheet != "lab") 
        assign(sheet, get(sheet) %>% replace_with_na_all(condition = ~.x %in% na_strings))
}

#
# Core variables all patients have data in
#
cv <- dplyr::select(corevar, all_of(c("SUBJID", "AGE", "SEX", "RACE", "PRHNTRTC", "B_ECOGCT", "DIAGTYCD", "TUMCAT", "PMAB", "HPVCD")))

# most the variable entries for EGFR related info is missing, no need to include
# for (colid in 1:dim(demo)[2]) {
#     print(colnames(demo)[colid])
#     print(length(which(is.na(demo[,colid]))) / dim(demo)[1])
# }

#
# Demographic variables
#
dem <- dplyr::select(demo, all_of(c("SUBJID", "B_WEIGHT", "B_HEIGHT", "B_BSA", "PRSURG", "PRRADIO", "DIAGMONS", "METMONS", "DIAGSDCD", "HDIFFMCD", "HSSBTMCD", "DISTMET", "CHILDPOT")))

#
# Baseline lab results
#
bl_lab <- filter(lab, LBBASE == "Y")
bl_lab <- bl_lab %>% replace_with_na_all(condition = ~.x %in% na_strings)
bl_lab <- dplyr::select(bl_lab, all_of(c("SUBJID", "LBTEST", "LBBLRES")))
bl_lab <- spread(bl_lab, LBTEST, LBBLRES)
colnames(bl_lab) <- gsub(" ", "_", colnames(bl_lab))

#
# Eligible patients
#
elig <- dplyr::select(eligcrit, all_of(c("SUBJID", "ELIGIBLE")))

#
# Baseline lesion
#
bl_ls <- filter(lesion, VISIT == "Screening")
bl_ls <- bl_ls %>% group_by(SUBJID) %>% group_by(SUBJID, ENRREFDY) %>% add_count(LSCAT, name = "LSCAT_count")  %>% group_by(SUBJID, LSCAT) %>% slice_max(ENRREFDY, n=1, with_ties = FALSE) %>% group_by(SUBJID, LSCAT) %>% mutate(avg_count = mean(LSCAT_count, na.rm = TRUE), avg_LSSLD = mean(LSSLD, na.rm = TRUE))

bl_ls <- bl_ls %>% group_by(SUBJID, LSCAT) %>% dplyr::slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(avg_count, avg_LSSLD, LSLD)) %>% select(-all_of(c("avg_LSSLD_Non-target lesion", "LSLD_Non-target lesion")))
colnames(bl_ls) <- c("SUBJID", "non_target_count", "target_count", "target_LSSLD", "target_LSLD")

#
# Medical history (only number of unresolved body system entries included)
#
# current_mh <- filter(medhist, MHSTAT == "Unresolved")
# current_mh <- current_mh %>% group_by(SUBJID) %>% add_count(MHBODSYS, name = "bodsys_count")
# current_mh <- dplyr::select(current_mh, all_of(c("SUBJID", "MHBODSYS", "bodsys_count")))
# current_mh <- unique(current_mh)
# current_mh <- spread(current_mh, MHBODSYS, bodsys_count) %>% dplyr::select(-"<NA>")
# current_mh[is.na(current_mh)] <- 0
# colnames(current_mh) <- gsub(" ", "_", colnames(current_mh))
# colnames(current_mh)[2:dim(current_mh)[2]] <- paste0("MH_", colnames(current_mh)[2:dim(current_mh)[2]])

#
# Other prior cancer history, radiotherapy or chemo therapy histories
#
chany <- dplyr::select(chemotx, all_of(c("SUBJID", "CHANY")))
chany <- unique(chany)
# other cancer tx history
otany <- dplyr::select(othcantx, all_of(c("SUBJID", "OTANY")))
otany <- unique(otany)

# radiotherapy history
# raany <- radiotx %>% group_by(SUBJID) %>% add_count(RASITE, name = "rasite_count")
# raany <- dplyr::select(raany, all_of(c("SUBJID", "RASITE", "rasite_count")))
# raany <- unique(raany)
# raany <- spread(raany, RASITE, rasite_count) %>% dplyr::select(-"<NA>")
# raany[is.na(raany)] <- 0
# colnames(raany) <- gsub(" ", "_", colnames(raany))
# colnames(raany)[2:dim(raany)[2]] <- paste0("RA_", colnames(raany)[2:dim(raany)[2]])

# Surgery history
surgh <- surghist %>% group_by(SUBJID) %>% add_count(SXTYPE, name = "surgtype_count")
surgh <- dplyr::select(surgh, all_of(c("SUBJID", "SXTYPE", "surgtype_count")))
surgh <- unique(surgh)
surgh <- spread(surgh, SXTYPE, surgtype_count) %>% dplyr::select(-"<NA>")


#
# Join and combine all baseline characteristics
#
cleaned_dat <- c("elig", "bl_ls", "bl_lab", "chany", "otany", "surgh")

X <- left_join(cv, dem, by = c("SUBJID"))
for (datf in cleaned_dat) {
    X <- left_join(X, get(datf), by = c("SUBJID"))
}
X <- unique(X) # some rows are duplicated after joining
identical(X$SUBJID, aeendpt$SUBJID) # check if the indexing is matched
X$SUBJID <- NULL
X <- missing_too_much(X, missing_threshold = 0.4)
X_imp <- impute_df_missing(clin_df = as.data.frame(X), save_ddt = FALSE)


identical(corevar$SUBJID, aeendpt$SUBJID) # check if the indexing is matched
W <- as.numeric(corevar$ATRT == "panit. plus chemotherapy")

X_imp <- as.data.frame(X_imp)
X_imp$ClNa <- X_imp$Chloride / X_imp$Sodium
X_imp$PWR <- X_imp$Platelets / X_imp$White_Blood_Cells
X_imp$ASTALT <- X_imp$Aspartate_Amino_Transferase / X_imp$Alanine_Amino_Transferase
X_imp$ASTLDH <- X_imp$Aspartate_Amino_Transferase / X_imp$Lactate_Dehydrogenase
X_imp$AGR <- X_imp$Albumin / X_imp$Hemoglobin
X_imp$NLR <- X_imp$Absolute_Neutrophil_Count / X_imp$Lymphocytes


vif_rmv <- c("B_HEIGHT", "B_BSA", "White_Blood_Cells", "Lymphocytes", "Total_Neutrophils", "Absolute_Neutrophil_Count", "B_WEIGHT", "Sodium", "Chloride", "Aspartate_Amino_Transferase", "Alanine_Amino_Transferase","Lactate_Dehydrogenase", "Red_Blood_Cells", "Hematocrit", "PRRADIO", "Albumin", "Hemoglobin", "Chloride", "Sodium")

X_imp <- dplyr::select(X_imp, -all_of(c(vif_rmv)))
X_imp <- as.matrix(X_imp)
#
# Export trial data
#
NCT00460265 <- list(X_imp, W)

# outcome measurements
OS <- data.frame(T = aeendpt$DTHDY / 30.4167, C = aeendpt$DTH) # primary outcome
PFS <- data.frame(T = aeendpt$PFSDYLRA / 30.4167, C = aeendpt$PFSLR) # secondary outcome
PFS[which(is.na(PFS$T)),1] <- OS[which(is.na(PFS$T)), 1]
RSP <- aeendpt$ROSLR # secondary outcome, changing from time until an event: CR or PR, is relative less informative than

# add AE
toxicity <- ae %>% group_by(SUBJID) %>% mutate(total_ae_grade = sum(AESEVCD)) %>% select(all_of(c("SUBJID", "total_ae_grade"))) %>% unique()
cv <- left_join(cv, toxicity, by = c("SUBJID"))
AE <- cv$total_ae_grade
AE[which(is.na(AE))] <- 0

NCT00460265_outcomes <- c("OS", "PFS", "RSP", "AE")

for (outcome in NCT00460265_outcomes) {
    if (outcome != "RSP" & outcome != "AE") {
    assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = get(outcome)[,2], X = X_imp)))
    } else {
        assign(paste0(outcome, "_Y_list"), list(get(outcome), NA, NA))
    }
}


save(NCT00460265, NCT00460265_outcomes, OS_Y_list, PFS_Y_list, RSP_Y_list, AE_Y_list, file = "./bin/load_RCT/RCT_obj/NCT00460265.RData")