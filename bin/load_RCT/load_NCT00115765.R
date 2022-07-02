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

datals <- c("corevar", "a_eendpt", "ae", "respeval","demo", "eligcrit", "lab", "lesion", "medhist", "chemotx", "radiotx", "surghist", "vitals_v")

for (sheet in datals) {
    assign(sheet, read.csv(paste0(file_path, sheet, ".csv")))
    if (sheet != "lab" & sheet != "lesion" & sheet != "medhist" & sheet != "vitals_v" & sheet != "ae")
        assign(sheet, get(sheet) %>% replace_with_na_all(condition = ~.x %in% na_strings))
}

#
# Core variables all patients have data in
#
cv <- dplyr::select(corevar, all_of(c("SUBJID", "ACHEMM", "AGE", "SEXCD", "RACCAT", "B_LESNM", "B_LDHYN", "B_LESNMT", "B_LESNMN", "B_METANM", "DURSURG", "PRADJYN", "B_ECOG", "DIAGTYPE", "PDADEN", "PDPROMED", "KRAS")))
# remove B_LDHYN because medical biochem test results include the continuous version already

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

lab_ddt <- dplyr::select(lab, all_of(c("LBTEST", "LBSTUNIT", "LBPANEL")))
lab_ddt <- unique(lab_ddt)

#
# Eligible patients
#
elig <- dplyr::select(eligcrit, all_of(c("SUBJID", "ELIGIBLE")))

# no need to include lesion data again

# NO NEED TO INCLUDE MH BECAUSE ALL ENTRIES ARE EVENTUALLY REMOVED DUE TO MISSINGNESS
#
# Medical history (only number of unresolved body system entries included)
#
# current_mh <- filter(medhist, MHSTAT == "Current or continuing")
# current_mh <- current_mh %>% group_by(SUBJID) %>% add_count(MHTERM, name = "cond_count")
# current_mh <- dplyr::select(current_mh, all_of(c("SUBJID", "MHTERM", "cond_count")))
# current_mh <- unique(current_mh)
# current_mh <- spread(current_mh, MHTERM, cond_count)
# current_mh[is.na(current_mh)] <- 0
# colnames(current_mh) <- gsub(" ", "_", colnames(current_mh))
# colnames(current_mh)[2:dim(current_mh)[2]] <- paste0("MH_", colnames(current_mh)[2:dim(current_mh)[2]])

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
cleaned_dat <- c("bl_lab", "elig", "chany", "raany", "surgh", "vtl")

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

# address VIF issues
X_imp$BUNCre <- X_imp$Blood_Urea_Nitrogen / X_imp$Creatinine
X_imp$NLR <- X_imp$Absolute_Neutrophil_Count / X_imp$Lymphocytes
X_imp$ASTALT <- X_imp$Aspartate_Amino_Transferase / X_imp$Alanine_Amino_Transferase
X_imp$ASTLDH <- X_imp$Aspartate_Amino_Transferase / X_imp$Lactate_Dehydrogenase
X_imp$ClNa <- X_imp$Chloride / X_imp$Sodium
X_imp$AGR <- X_imp$Albumin / X_imp$Hemoglobin

# Vitals
X_imp$MBP <- X_imp$Diastolic_blood_pressure + 1/3 * (X_imp$Systolic_blood_pressure - X_imp$Diastolic_blood_pressure)
# X_imp$BMI <- X_imp$B_WEIGHT / (X_imp$B_HEIGHT / 100)^2


vif_rmv <- c("Blood_Urea_Nitrogen", "Creatinine", "Absolute_Neutrophil_Count", "Lymphocytes", "Aspartate_Amino_Transferase", "Alanine_Amino_Transferase", "Lactate_Dehydrogenase", "Total_Neutrophils_(pct)", "Lymphocytes_(%)", "Red_Blood_Cells", "Albumin", "Hemoglobin", "Chloride", "Sodium", "Diastolic_blood_pressure", "Systolic_blood_pressure", "B_WEIGHT", "PRADJYN")

X_imp <- dplyr::select(X_imp, -all_of(c(vif_rmv)))
X_imp <- as.matrix(X_imp)

#
# Export trial data
#
NCT00115765 <- list(X_imp, W)

# Endpoints
endpt <- a_eendpt[a_eendpt$SUBJID %in% sel_SUBJID,]
OS <- data.frame(T = endpt$DTHMT, C = endpt$DTH) # primary outcome
PFS <- data.frame(T = endpt$PFSMTLR, C = endpt$PFSCR) # secondary outcome
TF <- data.frame(T = endpt$TFAILMT, C = endpt$TFAIL)

# Ordinal categorical outcome
RSP <- endpt$ROSLRCA # secondary outcome
rsps <- c("CR", "PR", "SD", "PD", "UE", "ND")
i <- 1
for (rsp_cat in rsps) {
    RSP[RSP == rsp_cat] <- i
    i <- i+1
}

# add AE
toxicity <- ae %>% group_by(SUBJID) %>% mutate(total_ae_grade = sum(AESEVCD)) %>% select(all_of(c("SUBJID", "total_ae_grade"))) %>% unique()
cv <- left_join(cv, toxicity, by = c("SUBJID" = "SUBJID"))
AE <- cv$total_ae_grade
AE[which(is.na(AE))] <- 0

NCT00115765_outcomes <- c("OS", "PFS", "TF", "RSP", "AE")

for (outcome in NCT00115765_outcomes) {
    if (!(outcome == "AE" | outcome == "RSP")) {
        sur_imp_res <- do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp))
        print(paste0("Infinite value check: ", outcome))
        for (entry in sur_imp_res) {
            print(which(is.infinite(entry)))
        }
        assign(paste0(outcome, "_Y_list"), sur_imp_res)
    } else {
        assign(paste0(outcome, "_Y_list"), list(get(outcome), NA, NA))
    }
}

save(NCT00115765, NCT00115765_outcomes, OS_Y_list, PFS_Y_list, TF_Y_list, RSP_Y_list, AE_Y_list, file = "./bin/load_RCT/RCT_obj/NCT00115765.RData")