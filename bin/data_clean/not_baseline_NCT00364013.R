library(dplyr)
library(impute)
library(tidyr)

adsl <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/adsl_pds2019.csv")

# check if assigned trt is the same as actual trt
all(adsl$TRT==adsl$ATRT)
message(c("Which patient(s) had different assigned vs actual treatment: ", which(!(adsl$TRT==adsl$ATRT))))

# disgard assigned treatment
adsl$TRT <- NULL
adsl[adsl == ""] <- NA

# arm 1: Panitumumab + FOLFOX will be considered as W=1, arm 2: FOLFOX alone will be considered as W=0
adsl$ATRT <- as.numeric(adsl$ATRT == "Panitumumab + FOLFOX")
print(table(adsl$ATRT))

##
## Append additional data sets to the primary adsl
##

#
# number of new target leison (ADLS)
#
# Adding leison (adls) data
adls <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/adls_pds2019.csv")
adls[adls == "NA"|adls == ""] <- NA
adls <- adls %>% filter(LSREADER != "Oncologist")

baseline_ls <- adls %>% filter(VISITDY < 0) %>% group_by(SUBJID) %>% group_by(SUBJID, VISITDY, LSREADER) %>% add_count(LSCAT, name = "LSCAT_count")  %>% group_by(SUBJID, LSCAT, LSREADER) %>% slice_max(VISITDY, n=1, with_ties = FALSE) %>% group_by(SUBJID, LSCAT) %>% mutate(avg_count = mean(LSCAT_count, na.rm = TRUE), avg_LSSLD = mean(LSSLD, na.rm = TRUE))

baseline_ls <- baseline_ls %>% group_by(SUBJID, LSCAT) %>% slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(avg_count, avg_LSSLD, LSLD)) %>% select(-all_of(c("avg_LSSLD_Non-target lesion", "LSLD_Non-target lesion")))
colnames(baseline_ls) <- c("SUBJID", "non_target_count", "target_count", "target_LSSLD", "target_LSLD")
# no-target lesions do not have diameters measured, so all LSLD values for non-target lesions entries are NAs

terminal_ls <- adls %>% filter(VISITDY > 0) %>% group_by(SUBJID) %>% group_by(SUBJID, VISITDY, LSREADER) %>% add_count(LSCAT, name = "LSCAT_count")  %>% group_by(SUBJID, LSCAT, LSREADER) %>% slice_max(VISITDY, n=1, with_ties = FALSE) %>% group_by(SUBJID, LSCAT) %>% mutate(avg_count = mean(LSCAT_count, na.rm = TRUE), avg_LSSLD = mean(LSSLD, na.rm = TRUE)) #LSCAT_count is the no of LS for this particular category of LS
terminal_ls <- terminal_ls %>% group_by(SUBJID, LSCAT) %>% slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(avg_count, avg_LSSLD, LSLD)) %>% mutate(new_lesion_count = replace(`avg_count_New lesion`, is.na(`avg_count_New lesion`), 0)) %>% select(-all_of(c("avg_count_New lesion", "avg_LSSLD_Non-target lesion", "avg_LSSLD_New lesion", "LSLD_Non-target lesion")))

colnames(terminal_ls) <- c("SUBJID", "non_target_count", "target_count", "target_LSSLD", "target_LSLD", "new_LSLD", "new_lesion_count")
# previously missing 
# sub_missing <- terminal_ls[apply(terminal_ls[,c(2:6)], 1, function(X) all(is.na(X))),1] # which patients do not have any terminal LS etnries

# Replace LSLD of new lesions with 0 if the patient has no new lesion detected
terminal_ls$new_LSLD[terminal_ls$new_lesion_count == 0] <- 0


# Separate lesion assements by oncologists/radiologists
#

oncologist_count <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/adls_pds2019.csv")
baseline_ls_onco <- oncologist_count %>% filter(LSREADER == "Oncologist", VISITDY < 0) %>% group_by(SUBJID, LSCAT) %>% slice_max(VISITDY, with_ties = TRUE) %>% add_count(LSCAT, name = "LSCAT_count") %>% group_by(SUBJID, LSCAT)  %>% slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(LSCAT_count, LSLD))
colnames(baseline_ls_onco) <- c("SUBJID", "non_target_count", "non_target_LSLD")

terminal_ls_onco <- oncologist_count %>% filter(LSREADER == "Oncologist", VISITDY > 0) %>% group_by(SUBJID, LSCAT) %>% slice_max(VISITDY, with_ties = TRUE) %>% add_count(LSCAT, name = "LSCAT_count") %>% group_by(SUBJID, LSCAT)  %>% slice(n=1) %>% pivot_wider(id_cols = SUBJID, names_from = LSCAT, values_from = c(LSCAT_count, LSLD)) %>% select(-all_of(c("LSLD_Non-target lesion")))

colnames(terminal_ls_onco) <- c("SUBJID", "new_lesion_count", "non_target_count", "new_LSLD")

# replace radiologist non-target counts with oncologists' non-target count
replace_subjid <- baseline_ls$SUBJID[baseline_ls$SUBJID %in% baseline_ls_onco$SUBJID]

for(id in replace_subjid) {
    onco_value <- baseline_ls_onco$non_target_count[baseline_ls_onco$SUBJID == id]
    radio_value <- baseline_ls$non_target_count[baseline_ls$SUBJID == id]

    if(onco_value > radio_value) {
        baseline_ls$non_target_count[baseline_ls$SUBJID == id] <- onco_value
    }
}

replace_subjid_terminal <- terminal_ls$SUBJID[terminal_ls$SUBJID %in% terminal_ls_onco$SUBJID]

for(id in replace_subjid_terminal) {
    onco_value <- terminal_ls_onco$non_target_count[terminal_ls_onco$SUBJID == id]
    radio_value <- terminal_ls$non_target_count[terminal_ls$SUBJID == id]

    if(is.na(radio_value) && is.na(onco_value))
        next
    else if(is.na(radio_value) | is.na(onco_value))
        terminal_ls$non_target_count[terminal_ls$SUBJID == id] <- ifelse(is.na(onco_value), radio_value, onco_value)
    else if(onco_value > radio_value)
        terminal_ls$non_target_count[terminal_ls$SUBJID == id] <- onco_value

    onco_nls <- terminal_ls_onco$new_lesion_count[terminal_ls_onco$SUBJID == id]
    radio_nls <- terminal_ls$new_lesion_count[terminal_ls$SUBJID == id]
    if(is.na(radio_nls))
        next
    else if(is.na(radio_nls) | is.na(onco_nls))
        terminal_ls$new_lesion_count[terminal_ls$SUBJID == id] <- ifelse(is.na(onco_nls), radio_nls, onco_nls)
    else 
        if(onco_nls > radio_nls)
        terminal_ls$new_lesion_count[terminal_ls$SUBJID == id] <- onco_nls

}

write.csv(baseline_ls, "./dat/NCT00364013_blls.csv", row.names = FALSE)
write.csv(terminal_ls, "./dat/NCT00364013_tlls.csv", row.names = FALSE)

#
# Expert rated response (ADRSP)
#
adrsp <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/adrsp_pds2019.csv")


# during screening, no RSRESP & no new RSNEW were recorded
# respnose status by final visit day
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

write.csv(rsp, "./dat/NCT00364013_rsp.csv", row.names = FALSE)

# cross validate new lesions data with adls
new_ls_sub <- unique(adrsp$SUBJID[adrsp$RSNEW == "Y"])
check_ls <- terminal_ls %>% filter(new_lesion_count > 0)
length(new_ls_sub)
dim(check_ls)
message(c("Did ADLS capture all new lesions identifiable by experts? ", all(new_ls_sub %in% check_ls$SUBJID)))
missing_ls <- new_ls_sub[!(new_ls_sub %in% check_ls$SUBJID)]
rsp[rsp$SUBJID %in% missing_ls,]


#
# Adverse effect reports (ADAE)
#
adae <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/adae_pds2019.csv")

# patients with continuing adverse events will have their AEENDYI be replaced by their last follow-up day in the DTHDY column of adsl
cont_sub <- data.frame(cont_subjid = adae$SUBJID[adae$AECONT == "Y"])
cont_sub <- left_join(cont_sub, adsl, by = c("cont_subjid" = "SUBJID"))
cont_sub <- cont_sub %>% select(all_of(c("cont_subjid", "DTHDY"))) %>% rename(SUBJID = cont_subjid)

sum_ae <- left_join(adae, cont_sub, by = "SUBJID")

sum_ae <- sum_ae %>% 
            mutate(AEENDYI = coalesce(AEENDYI, DTHDY)) %>% 
            group_by(SUBJID) %>% 
            mutate(AECONT = ifelse(AECONT == "Y", 1, 0)) %>% 
            mutate(ae_duration_sum = ifelse(AECONT == 0, sum(AEENDYI - AESTDYI), sum(DTHDY - AESTDYI))) %>% 
            add_tally(AESEVCD, name = "ae_severity_count") %>% 
            add_tally(AECONT, name = "ae_contd_count") %>% 
            select(-all_of(c("AEPT", "AESOC", "AESTDYI", "AEENDYI", "AESEVCD", "AECONT", "DTHDY"))) %>% 
            slice(n = 1)
write.csv(sum_ae, "./dat/NCT00364013_ae_sum.csv", row.names = FALSE)

#
# Biomarker statuses
#
biomark <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_309_NCT00364013/csv/biomark_pds2019.csv")

biomarknm <- unlist(biomark[1,seq(2, dim(biomark)[2], 2)])
biomarknm <- gsub(" ", "_", biomarknm)

biom <- biomark[,c(1,seq(3, dim(biomark)[2], 2))]
colnames(biom) <- c("SUBJID", biomarknm)
biom[biom == ""] <- NA

write.csv(biom, "./dat/NCT00364013_biom.csv", row.names = FALSE)

#######################################
#######################################
#######################################
clin_df <- adsl
# convert categorical to numeric
num_df <- cat2num(clin_df)
num_clin_df <- num_df[[1]]
clin_df_ddt <- num_df[[2]]

#
# impute missing data if necessary
#
if (any(is.na(num_clin_df))) {
    pre_proc <- scale(num_clin_df, center = TRUE, scale = TRUE)
    imp_clin_df <- impute.knn(pre_proc, k=10)
    message(c("Seed used: ", imp_clin_df$rng.seed))
    imp_clin_df <- t(apply(imp_clin_df$data, 1, function(r) r *attr(imp_clin_df$data,'scaled:scale') + attr(imp_clin_df$data, 'scaled:center')))
}
message("Is there still missing data? ", any(is.na(imp_clin_df)))
imp_clin_df <- as.data.frame(imp_clin_df)
class(imp_clin_df$PFSCR) <- "integer" # after imputation returned floating point for 0s

write.csv(imp_clin_df, "./dat/NCT00364013_adls.csv", row.names = FALSE)

# clin_df$SUBJID <- NULL # remove subject id, note that only 80% of the patients were included in the final dataset