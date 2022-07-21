#' Load dataset for NCT00364013
#' Primary outcome: Progressive disease defined as least a 20% increase in the sum of the longest diameters (SLD) of target lesions, taking as reference the nadir SLD recorded since the treatment started or the appearance of one or more new lesions, or the unequivocal progression of existing non-target lesions. 
#' Secondary outcome: 
#' - Overall Survival
#' -  Objective Tumor Response
#' - Time to Treatment Failure 
#' - Progression-free Survival Time (Wild-type KRAS)
#' -  Progression-free Survival Time (Mutant KRAS)

file_path <- "./dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/"

datals <- c("adsl", "adlb", "adls", "adrsp", "biomark", "adae")

# set to NA in R language for any type of missing or not otherwise specified entries
for (sheet in datals) {
    assign(sheet, read.csv(paste0(file_path, sheet, "_pds2019.csv")))
    if (sheet != "adlb")
        assign(sheet, get(sheet) %>% replace_with_na_all(condition = ~.x %in% na_strings))
}

adsl$SUBJID[which(is.na(adsl$SUBJID))] <- c(98, 99)
corevar <- dplyr::select(adsl, -all_of(c("TRT", "ATRT", "DTHDY", "DTH", "PFSDYCR", "PFSCR")))

#
# Baseline lab results
#
bl_lab <- filter(adlb, LBBASE == "Y") %>% group_by(SUBJID) %>%  slice_max(VISITDY, n=1, with_ties = TRUE)
bl_lab <- dplyr::select(bl_lab, all_of(c("SUBJID", "LBTEST", "LBSTRESN")))
bl_lab <- spread(bl_lab, LBTEST, LBSTRESN)
colnames(bl_lab) <- gsub(" ", "_", colnames(bl_lab))
bl_lab <- bl_lab %>% replace_with_na_all(condition = ~.x %in% na_strings)
# remove entries failed to be identified from SUBJID
bl_lab <- bl_lab[which(!is.na(bl_lab$SUBJID)),]

#
# Baseline lesion
#
bl_ls <- filter(adls, VISIT == "Screening")
bl_ls <- bl_ls[which(!is.na(bl_ls$SUBJID)),]
bl_ls <- bl_ls %>% group_by(SUBJID) %>% group_by(SUBJID, VISITDY) %>% add_count(LSCAT, name = "LSCAT_count")  %>% group_by(SUBJID, LSCAT) %>% slice_max(VISITDY, n=1, with_ties = FALSE) %>% group_by(SUBJID, LSCAT) %>% mutate(avg_count = mean(LSCAT_count, na.rm = TRUE), avg_LSSLD = mean(LSSLD, na.rm = TRUE))

bl_ls <- bl_ls %>% group_by(SUBJID, LSCAT) %>% dplyr::slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(avg_count, avg_LSSLD, LSLD)) %>% select(-all_of(c("avg_LSSLD_Non-target lesion", "LSLD_Non-target lesion")))
colnames(bl_ls) <- c("SUBJID", "non_target_count", "target_count", "target_LSSLD", "target_LSLD")

#
# Baseline biomarker
#
biomarknm <- unlist(biomark[1,seq(2, dim(biomark)[2], 2)])
biomarknm <- gsub(" ", "_", biomarknm)

biom <- biomark[,c(1,seq(3, dim(biomark)[2], 2))]
colnames(biom) <- c("SUBJID", biomarknm)
biom$SUBJID[which(is.na(biom$SUBJID))] <- c(98, 99)
biom <- biom %>% group_by(SUBJID) %>% fill(everything(), .direction = "downup") %>% dplyr::slice(1)
# biom[is.na(biom)] <- "Unknown"

# only for checking missing proportion
# for (col in colnames(biom)) {
#     print(length(which(is.na(biom[[col]]))) / dim(biom)[1])
# }

#
# Join and combine all baseline characteristics
# 
cleaned_dat <- c("bl_lab","bl_ls")

X <- left_join(corevar, biom, by = c("SUBJID"))
for (datf in cleaned_dat) {
    X <- left_join(X, get(datf), by = c("SUBJID"))
    print(dim(unique(X)))
}
sel_SUBJID <- X$SUBJID
X$SUBJID <- NULL
X <- missing_too_much(X, missing_threshold = 0.5)
X_imp <- impute_df_missing(clin_df = as.data.frame(X), save_ddt = FALSE)

# check if assigned trt is the same as actual trt
all(adsl$TRT==adsl$ATRT)
message(c("Which patient(s) had different assigned vs actual treatment: ", which(!(adsl$TRT==adsl$ATRT))))

W <- adsl$ATRT[adsl$SUBJID %in% sel_SUBJID]
W <- as.numeric(W == "Panitumumab + FOLFOX")

# Transformation

X_imp$AGR <- X_imp$Albumin / X_imp$Hemoglobin
X_imp$PWR <- X_imp$Platelets / X_imp$White_Blood_Cells

vif_rmv <- c("Albumin", "Hemoglobin", "Platelets", "White_Blood_Cells")

X_imp <- dplyr::select(X_imp, -all_of(c(vif_rmv)))
X_imp <- as.matrix(X_imp)

NCT00364013 <- list(X_imp, W)

# outcome measurements
adsl <- adsl[adsl$SUBJID %in% sel_SUBJID,]
outcomes <- dplyr::select(adsl, all_of(c("SUBJID", "DTHDY", "DTH", "PFSDYCR", "PFSCR")))

# Objective response
terminal_rsp <- adrsp %>% group_by(SUBJID) %>% slice_max(VISITDY, n=1, with_ties = TRUE) %>% pivot_wider(id_cols = SUBJID, names_from = RSREADER, values_from = RSRESP)
terminal_rsp <- as.data.frame(terminal_rsp)
# encode responses in graded categorical format
rsps <- c("Complete response", "Partial response", "Stable disease", "Progressive disease")
i <- 1
for (rsp_cat in rsps) {
    terminal_rsp[terminal_rsp == rsp_cat] <- i
    i <- i+1
}
rsp <- as.data.frame(apply(terminal_rsp, 2, as.numeric))
rsp_avg <- apply(rsp[,c(2:4)], 1, function(X) mean(X, na.rm = TRUE))
rsp <- cbind(rsp, rsp_avg)
rsp <- select(rsp, all_of(c("SUBJID", "rsp_avg")))

# Impute missing responses
imp_df <- left_join(outcomes, rsp, by = c("SUBJID"))
imp_outcome <- impute_df_missing(clin_df = as.data.frame(imp_df[,2:6]), save_ddt = FALSE)
imp_df <- cbind(data.frame(SUBJID = imp_df$SUBJID), imp_outcome)

OS <- data.frame(T = imp_df$DTHDY * 0.142857, C = imp_df$DTH) # primary outcome
PFS <- data.frame(T = imp_df$PFSDYCR * 0.142857, C = imp_df$PFSCR) # secondary outcome
RSP <- imp_df$rsp_avg

# add AE
toxicity <- adae %>% group_by(SUBJID) %>% mutate(total_ae_grade = sum(AESEVCD)) %>% select(all_of(c("SUBJID", "total_ae_grade"))) %>% unique()
corevar <- left_join(corevar, toxicity, by = c("SUBJID"))
AE <- corevar$total_ae_grade
AE[which(is.na(AE))] <- 0

NCT00364013_outcomes <- c("OS", "PFS", "RSP", "AE")
 

for (outcome in NCT00364013_outcomes) {
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

 
save(NCT00364013, NCT00364013_outcomes, OS_Y_list, PFS_Y_list, RSP_Y_list, AE_Y_list, file = "./bin/load_RCT/RCT_obj/NCT00364013.RData")