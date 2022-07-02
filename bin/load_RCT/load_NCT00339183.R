#' Script to load the trial NCT00339183
#' Pirmary objective: Progression-free Survival (PFS) and Overall Survival. Both stratified by KRAS mutation status
#' Secondary objective: 
#' 
#' 

file_path <- "./dat/PDS/Colorec_Amgen_2006_263_NCT00339183/csv/"

datals <- c("corevar", "a_eendpt", "a_endpt", "ae", "demo", "eligcrit", "lab", "lesion", "chemotx", "othcantx", "surghist", "vitals_v", "disposit")

# set to NA in R language for any type of missing or not otherwise specified entries
# na strings loaded in the load_lib.R script
for (sheet in datals) {
    assign(sheet, read.csv(paste0(file_path, sheet, ".csv")))
    if (sheet != "lab" & sheet != "lesion" & sheet != "medhist" & sheet != "vitals_v") 
        assign(sheet, get(sheet) %>% replace_with_na_all(condition = ~.x %in% na_strings))
}

#
# Core variables all patients have data in
#
cv <- dplyr::select(corevar, all_of(c("SUBJID", "AGE", "SEXCD", "RACCAT", "PRBEVCDR", "PRBEVCDM", "PROXACDR", "PROXACDM", "B_METANM", "LIVRONLY", "PRADJYN", "DPADJUV", "B_LDHYN", "B_LDH2YN", "B_LDHNM", "B_ECOGI", "DIAGTYPE", "KRAS")))

#
# Demographic variables
#
dem <- dplyr::select(demo, all_of(c("SUBJID", "B_WEIGHT", "B_HEIGHT", "B_BSA", "DIAGMONS", "METMONS", "HISTYPE", "HISSUBTY", "HDIFFER", "BIOPTPCT", "BIOPRS", "BIOPS123", "BIOPMAX", "BIOPPCT", "BIOPCYTO", "BIOPMAXC", "CHILDPOT")))

# for (col in colnames(demo)) {
#     print(col)
#     print(length(which(is.na(demo[[col]]))) / dim(demo)[1])
# }

#
# Eligible patients
#
elig <- dplyr::select(eligcrit, all_of(c("SUBJID", "ELIGIBLE")))

#
# Baseline lab results
#
bl_lab <- filter(lab, LBBASE == "Y")
bl_lab <- dplyr::select(bl_lab, all_of(c("SUBJID", "LBTEST", "LBBLRES")))
bl_lab <- bl_lab %>% replace_with_na_all(condition = ~.x %in% na_strings)
bl_lab <- spread(bl_lab, LBTEST, LBBLRES)
colnames(bl_lab) <- gsub(" ", "_", colnames(bl_lab))


#
# Baseline lesion
#
bl_ls <- filter(lesion, VISIT == "Screening")
bl_ls <- bl_ls %>% replace_with_na_all(condition = ~.x %in% na_strings)
bl_ls <- bl_ls %>% group_by(SUBJID) %>% group_by(SUBJID, ENRREFDY) %>% add_count(LSCAT, name = "LSCAT_count")  %>% group_by(SUBJID, LSCAT) %>% slice_max(ENRREFDY, n=1, with_ties = FALSE) %>% group_by(SUBJID, LSCAT) %>% mutate(avg_count = mean(LSCAT_count, na.rm = TRUE), avg_LSSLD = mean(LSSLD, na.rm = TRUE))

bl_ls <- bl_ls %>% group_by(SUBJID, LSCAT) %>% dplyr::slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(avg_count, avg_LSSLD, LSLD)) %>% select(-all_of(c("avg_LSSLD_Non-target lesion", "LSLD_Non-target lesion", "avg_count_NA", "avg_LSSLD_NA", "LSLD_NA")))
colnames(bl_ls) <- c("SUBJID", "non_target_count", "target_count", "target_LSSLD", "target_LSLD")


# Medical history (only number of unresolved body system entries included)
#
# current_mh <- filter(medhist, MHSTAT == "Current or continuing") %>% replace_with_na_all(condition = ~.x %in% na_strings)
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
cleaned_dat <- c("bl_lab", "bl_ls", "elig", "chany", "otany", "surgh", "vtl")

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

W <- as.numeric(corevar$ATRT == "Panitumumab + FOLFIRI")

# address VIF issues
# More than 30% missing data! Removing from df.: Blood_Urea_Nitrogen"

X_imp$NLR <- X_imp$Absolute_Neutrophil_Count / X_imp$Lymphocytes
X_imp$ASTALT <- X_imp$Aspartate_Amino_Transferase / X_imp$Alanine_Amino_Transferase
X_imp$ASTLDH <- X_imp$Aspartate_Amino_Transferase / X_imp$B_LDHNM
X_imp$ClNa <- X_imp$Chloride / X_imp$Sodium
X_imp$AGR <- X_imp$Albumin / X_imp$Hemoglobin
X_imp$PWR <- X_imp$Platelets / X_imp$Lymphocytes

# Vitals
X_imp$MBP <- X_imp$Diastolic_blood_pressure + 1/3 * (X_imp$Systolic_blood_pressure - X_imp$Diastolic_blood_pressure)
# X_imp$BMI <- X_imp$B_WEIGHT / (X_imp$B_HEIGHT / 100)^2


vif_rmv <- c("Absolute_Neutrophil_Count", "Lymphocytes", "Aspartate_Amino_Transferase", "Alanine_Amino_Transferase", "B_LDHNM", "B_LDHYN", "B_LDH2YN", "Lactate_Dehydrogenase", "Total_Neutrophils", "Lymphocytes", "Red_Blood_Cells", "Albumin", "Hemoglobin", "Chloride", "Sodium", "Diastolic_blood_pressure", "Systolic_blood_pressure", "B_WEIGHT", "Platelets")

X_imp <- dplyr::select(X_imp, -all_of(c(vif_rmv)))
X_imp <- as.matrix(X_imp)


#
# Export trial data
#
NCT00339183 <- list(X_imp, W)

# Objective response
endpt <- a_eendpt[a_eendpt$SUBJID %in% sel_SUBJID,]
endpt <- left_join(endpt, disposit, by = c("SUBJID"))
PFS <- data.frame(T = endpt$PDDYLR / 30.4167, C = endpt$PDLR) # secondary outcome
OS <- data.frame(T = endpt$DTHDY / 30.4167 , C = endpt$DTH) # primary outcome

# one DTHDY entry was missing, last contact day relative to 1st dose will be used instead
OS[which(is.na(OS$T)), 1] <- endpt$PFUP[which(is.na(OS$T))] / 7 # convert days to months


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
cv <- left_join(cv, toxicity, by = c("SUBJID"))
AE <- cv$total_ae_grade
AE[which(is.na(AE))] <- 0

NCT00339183_outcomes <- c("OS", "PFS", "RSP", "AE")

for (outcome in NCT00339183_outcomes) {
    if (outcome != "RSP" & outcome != "AE") {
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

save(NCT00339183, NCT00339183_outcomes, OS_Y_list, PFS_Y_list, RSP_Y_list, AE_Y_list, file = "./bin/load_RCT/RCT_obj/NCT00339183.RData")