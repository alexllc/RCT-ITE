
#' Primary outcome: Progression-free Survival Time [ Time Frame: From randomization to the data cut-off date of 30 June 2005. The median follow-up time was 20.0 weeks in the panitumumab plus BSC group and 18.2 weeks in the BSC alone group. ]
#' Secondary outcome: OS, Objective Tumor Response, Time to Treatment Failure, Duration of Stable Disease
file_path <- "./dat/PDS/Colorec_Amgen_2004_310_NCT00113763/csv/"

datals <- c("adsl", "adlb", "adls", "adrsp", "biomark")

# set to NA in R language for any type of missing or not otherwise specified entries
for (sheet in datals) {
    assign(sheet, read.csv(paste0(file_path, sheet, "_pds2019.csv")))
    if (sheet != "adlb")
        assign(sheet, get(sheet) %>% replace_with_na_all(condition = ~.x %in% na_strings))
}
# between SUBJID 97 to 101, there are two subjects with NA, they were simply named 98 and 99 but there were no way to identify the corresponding records from the lesion or the lab test dataset for these two
adsl$SUBJID[which(is.na(adsl$SUBJID))] <- c(98, 99)
corevar <- dplyr::select(adsl, -all_of(c("TRT", "ATRT", "DTHDYX", "DTHX", "PFSDYCR", "PFSCR")))

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

bl_ls <- bl_ls %>% group_by(SUBJID, LSCAT) %>% slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(avg_count, avg_LSSLD, LSLD)) %>% select(-all_of(c("avg_LSSLD_Non-target lesion", "LSLD_Non-target lesion")))
colnames(bl_ls) <- c("SUBJID", "non_target_count", "target_count", "target_LSSLD", "target_LSLD")
bl_ls <- unique(bl_ls)


#
# Baseline biomarker
#
biomarknm <- unlist(biomark[1,seq(2, dim(biomark)[2], 2)])
biomarknm <- gsub(" ", "_", biomarknm)

biom <- biomark[,c(1,seq(3, dim(biomark)[2], 2))]
colnames(biom) <- c("SUBJID", biomarknm)
biom$SUBJID[which(is.na(biom$SUBJID))] <- c(98, 99)
biom <- biom %>% group_by(SUBJID) %>% fill(everything(), .direction = "downup") %>% slice(1)

#
# Join and combine all baseline characteristics
#
cleaned_dat <- c("biom", "bl_ls")

X <- left_join(corevar, bl_lab, by = c("SUBJID"))
for (datf in cleaned_dat) {
    X <- left_join(X, get(datf), by = c("SUBJID"))
    print(dim(unique(X)))
}
sel_SUBJID <- X$SUBJID
X$SUBJID <- NULL
X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = as.data.frame(X), save_ddt = FALSE)

# check if assigned trt is the same as actual trt
all(adsl$TRT==adsl$ATRT)
message(c("Which patient(s) had different assigned vs actual treatment: ", which(!(adsl$TRT==adsl$ATRT))))

W <- adsl$ATRT[adsl$SUBJID %in% sel_SUBJID]
#
# Export trial data
#
NCT00113763 <- list(X_imp, W)

# outcome measurements
adsl <- adsl[adsl$SUBJID %in% sel_SUBJID,]
outcomes <- dplyr::select(adsl, all_of(c("SUBJID", "DTHDYX", "DTHX", "PFSDYCR", "PFSCR")))

# Objective response
terminal_rsp <- adrsp %>% group_by(SUBJID) %>% slice_max(VISITDY, n=1, with_ties = TRUE) %>% pivot_wider(id_cols = SUBJID, names_from = RSREADER, values_from = RSRESP)
terminal_rsp[terminal_rsp == "Unable to evaluate"] <- NA

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

OS <- data.frame(T = imp_df$DTHDYX / 30.417, C = imp_df$DTHX) # primary outcome
PFS <- data.frame(T = imp_df$PFSDYCR / 30.417, C = imp_df$PFSCR) # secondary outcome

NCT00113763_outcomes <- c("OS", "PFS", "RSP")
 

for (outcome in NCT00113763_outcomes) {
    if (outcome != "RSP") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp)))
    } else {
        RSP_Y_list <- list(imp_df$rsp_avg, NA, NA)
    }
}
