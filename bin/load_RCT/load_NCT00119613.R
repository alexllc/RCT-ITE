source("./bin/load_lib.R")

file_path <- "./dat/PDS/LungSm_Amgen_2002_266_NCT00119613/csv/"

datals <- list.files(path = file_path)

for (file in datals) {
    print(paste0("Processing: ", file))
    datas <- read.csv(paste0(file_path, file))
    sheetname <- strsplit(file, "_|\\.")[[1]][2]
    assign(sheetname, datas)
    if (sheetname != "lab" & sheetname != "lesion" & sheetname != "medhist" & sheetname != "vitals") 
        assign(sheetname, get(sheetname) %>% replace_with_na_all(condition = ~.x %in% na_strings))
}


#
# Key variables
#
kv <- dplyr::select(keyvar, all_of(c("SUBJID", "AGE", "SEX", "RACE", "B_HGB", "B_SEREPO", "TUMORCD", "B_ECOGN", "B_LDH2")))
kv$B_HGB[which(is.na(kv$B_HGB))] <- eendpt$HGB2[which(is.na(kv$B_HGB))]
kv <- kv[which(!is.na(kv$SUBJID)),]

# Load enrollment data
elig <- dplyr::select(entry, all_of(c("SUBJID", "ELIGYN")))
elig <- elig[which(!is.na(elig$SUBJID)),]
elig <- unique(elig)

#
# body characterisitcs
#
bcha <- dplyr::select(bchar, all_of(c("SUBJID", "B_HEIGHT", "B_WEIGHT", "PRGPOTYN")))
bcha$PRGPOTYN <- as.numeric(gsub("7|8|9", 0, bcha$PRGPOTYN))

#
# Process chemo dataset
#
chemo_sum <- chemo %>% group_by(SUBJID, TOTDOSE) %>% slice_max(CYCLE, n = 1, with_ties = FALSE) %>% group_by(SUBJID) %>% mutate(sum_totdoses = sum(TOTDOSE, na.rm = TRUE)) %>% select(all_of(c("SUBJID", "sum_totdoses"))) %>% unique()

#
# medical history
#
# Medical history (only number of unresolved body system entries included)
#
current_mh <- filter(medhis, MHSTAT == "Continuing") %>% replace_with_na_all(condition = ~.x %in% na_strings)
current_mh <- current_mh %>% group_by(SUBJID) %>% add_count(MEDHX, name = "sys_count") %>% dplyr::select(all_of(c("SUBJID", "MEDHX", "sys_count"))) %>% unique() %>% spread(MEDHX, sys_count)

current_mh[is.na(current_mh)] <- 0
colnames(current_mh) <- gsub(" ", "_", colnames(current_mh))
colnames(current_mh)[2:dim(current_mh)[2]] <- paste0("MH_", colnames(current_mh)[2:dim(current_mh)[2]])


#
# Transfusion
#
tf <- trans %>% group_by(SUBJID) %>% mutate(sum_trans = sum(VOLBLD, na.rm = TRUE)) %>% select(all_of(c("SUBJID", "sum_trans"))) %>% unique()

#
# Total doses of conmed taken, but expect most counts to be 0s
#
cmed <- filter(conmed, CMPRIOR == 1) %>% group_by(SUBJID) %>% add_count(DCTERM, name = "drug_class_count") %>% dplyr::select(all_of(c("SUBJID", "DCTERM", "drug_class_count"))) %>% unique() %>% spread(DCTERM, drug_class_count)

cmed[is.na(cmed)] <- 0
colnames(cmed) <- gsub(" ", "_", colnames(cmed))
colnames(cmed)[2:dim(cmed)[2]] <- paste0("MH_", colnames(cmed)[2:dim(cmed)[2]])

#
# vitals
#
vtl <- filter(vitals, PHASE == "Pre-treatment")
vtl <- vtl %>% replace_with_na_all(condition = ~.x %in% na_strings)
vtl <- vtl %>% group_by(SUBJID) %>% slice_max(STUDYDAY, n=1, with_ties = FALSE) %>% dplyr::select(all_of(c("SBP", "DBP", "PULSE", "ECOGCD")))

#
# hematology
#
hema <- filter(hemat, PHASE == "Pre-treatment") %>% group_by(SUBJID) %>% slice_max(STUDYDAY, n=1, with_ties = FALSE) %>% dplyr::select(all_of(c("SUBJID", "RBC", "HGB", "HCT", "PLT", "WBC", "ANC", "NEUT", "BANDS", "EOSIN", "SEGS")))

#
# Medical biochemistry (NOW DISCARDED IN FAVOR OF THE LAB DATASET)
#
# biochem <- filter(chem, PHASE == "Pre-treatment") %>% group_by(SUBJID) %>% slice_max(STUDYDAY, n=1, with_ties = FALSE) %>% dplyr::select(all_of(c("SUBJID", "SODIUM", "POTAS", "CL", "BICARB", "TLPROT", "ALBUMN", "GLUC", "BUN", "UREA", "CREAT", "URACID", "TLBILI", "ALKPH", "LDH", "SGOT", "SGPT")))
# colnames(biochem)[2:dim(biochem)[2]] <- paste0("BCHEM_", colnames(biochem)[2:dim(biochem)[2]])


#
# Iron
#
fe2 <- iron %>% replace_with_na_all(condition = ~.x %in% na_strings) %>% group_by(SUBJID) %>% slice_max(STUDYDAY, n=1, with_ties = FALSE) %>% summarise_all(coalesce_by_column) %>% dplyr::select(all_of(c("SUBJID", "SERMFE", "SERMTF", "FERRIT", "TFSAT", "TIBC", "FOLATE", "VITB12")))

# No baseline information for lesion data

# Head body scan dataset useless

#
# Radio
#
radio_sum <- radio %>% filter(PHASE == "Pre-treatment") %>% group_by(SUBJID) %>% mutate(sum_totdoses = sum(TOTDOSE, na.rm = TRUE)) %>% select(all_of(c("SUBJID", "sum_totdoses"))) %>% unique()

#
# Baseline lab results
#
bl_lab <- filter(lab, PHASE == "Pre-treatment") %>% replace_with_na_all(condition = ~.x %in% na_strings)
bl_lab <- bl_lab[which(!is.na(bl_lab$SUBJID)),]
bl_lab <- bl_lab %>% group_by(SUBJID) %>% slice() %>% dplyr::select(all_of(c("SUBJID", "TEST", "RESULT_B"))) %>% unique()
bl_lab <- spread(bl_lab, TEST, RESULT_B)
colnames(bl_lab)[2:dim(bl_lab)[2]] <- paste0("LB_", colnames(bl_lab)[2:dim(bl_lab)[2]])

#
# Join and combine all baseline characteristics
#
cleaned_dat <- c("bcha", "chemo_sum", "current_mh", "tf", "cmed", "vtl", "hema", "fe2", "radio_sum", "bl_lab")

X <- left_join(kv, elig, by = c("SUBJID"))
for (datf in cleaned_dat) {
    X <- left_join(X, get(datf), by = c("SUBJID"))
    print(dim(X))
}
X <- unique(X) # some rows are duplicated after joining
SUBJID_sel <- X$SUBJID
X$SUBJID <- NULL
X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = as.data.frame(X), save_ddt = FALSE)

keyvar <- keyvar[keybar$SUBJID %in% SUBJID_sel,]
W <- as.numeric(keyvar$TXGROUP == "NESP")

# Export dataset
NCT00119613 <- list(X_imp, W)

# Outcomes

# Primary: Change in hemoglobin concentration from baseline to the end of the chemotherapy treatment period 
# Primary: (Overall) Survival time
# Secondary: change in FACT-fatigue score
# Secondary: incidence of adverse events
# changes in laboratory values 
outcomes <- dplyr::select(eendpt, all_of(c("SUBJID", "CHHB", "DTHDY", "DTH", "PFSCD", "PFSDY")))
outcomes <- left_join(outcomes, dplyr::select(qol, all_of(c("SUBJID", "FATCHG1"))), by = c("SUBJID"))
outcomes <- outcomes[outcomes$SUBJID %in% SUBJID_sel,]

imp_outcome <- impute_df_missing(clin_df = as.data.frame(outcomes[,2:dim(outcomes)[2]]), save_ddt = FALSE)
imp_df <- cbind(data.frame(SUBJID = outcomes$SUBJID), imp_outcome)

OS <- data.frame(T = imp_df$DTHDY / 7, C = imp_df$DTH) # primary outcome
PFS <- data.frame(T = imp_df$PFSDY / 7, C = imp_df$PFSCD) # secondary outcome
CHHB <- imp_df$CHHB
FATCHG <- imp_df$FATCHG1
NCT00119613_outcomes <- c("OS", "PFS", "CHHB", "FATCHG")
 

for (outcome in NCT00119613_outcomes) {
    if (outcome != "CHHB" & outcome != "FATCHG") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp)))
    } else {
        assign(paste0(outcome, "_Y_list"), list(imp_df[[outcome]], NA, NA))
    }
}
