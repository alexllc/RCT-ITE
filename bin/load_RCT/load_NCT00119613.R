source("./bin/load_lib.R")
source("./bin/impute_survival.R")

file_path <- "./dat/PDS/LungSm_Amgen_2002_266_NCT00119613/csv/"

mode <- "core" # "extended"


if (mode == "extended") {
    ds_list <- c("c_keyvar", "c_bchar", "c_medhis", , "c_entry", "c_vitals", "c_hemat", "c_chem", "c_iron", "c_lesion", "c_hbscan", "a_lab")

    # Too little pre-treatment entries:
    # "c_radio", "c_chemo","c_conmed"

    for (file in ds_list) {
        imp <- read.csv(paste0(file_path, file, ".csv"))
        assign(strsplit(file, "_")[[1]][2], imp)
    }

    for (file in ds_list) {
        print(paste0(rep("=", 80), collapse = ""))
        file_name <- strsplit(file, "_")[[1]][2]
        print(paste0("Evaluating: ", file_name))

        count_entry <- filter(get(file_name), "STUDYDAY" < 0)
        print(table(get(file_name)[["STUDYDAY"]]))
    }

    # select extra variables in hemat, chem, lesion, hbscan

    var_sel <- list(c("AGE", "SEXCD", "RACECD", "B_HGB", "B_SEREPO", "TUMORCD", "B_ECOGN", "B_LDHN"), # keybar
                    c("B_WEIGHT", "B_HEIGHT"), # bchar
                    c("MEDHXCD", "MHSTATCD"), # medhis
                    c("RBC", "HGB", "HCT", "PLT", "WBC", "ANC", "NEUT", "BANDS", "EOSIN", "SEGS"), #hemat
                    c("SODIUM", "POTAS", "CL", "BICARB", "TLPROT", "ALBUMN", "GLUC", "BUN", "UREA", "CREAT", "URACID", "TLBILI", "ALKPH", "LDH", "SGOT", "SGPT"), # chem
                    c("CTSINDEX", "NTLALL"), # lesion
                    c("BONERCD", "HEADRCD")) # hbscan

    sel_list <- c("keyvar", "bchar", "medhis", "hemat", "chem", "lesion", "hbscan")

    for (i in 1:length(sel_list)) {
        sel_df <- get(sel_list[i])
        pre_treat_df <- try(filter(sel_df, STUDYDAY < 0))
        if (class(pre_treat_df) == "try-error") {
            print(paste0(sel_list[i], " Dataset does not have STUDYDAY"))
            sel_df <- dplyr::select(sel_df, all_of(c("SUBJID", var_sel[[i]])))
        } else if ( dim(pre_treat_df)[1] == 0) {
            sel_df <- dplyr::select(sel_df, all_of(c("SUBJID", var_sel[[i]])))
        } else {
            sel_df <- sel_df %>% group_by(SUBJID) %>% slice_max(STUDYDAY, n = 1)
        }
        assign(paste0("sel_", sel_list[i]), sel_df)
    }

    for (i in 1:length(sel_list)) {
        if (i == 1) {
            X <- left_join(get(paste0("sel_", sel_list[1])), get(paste0("sel_", sel_list[2])), by = c("SUBJID"))
        } else if (i == 2) {
            next
        } else {
            X <- left_join(X, get(paste0("sel_", sel_list[i])), by = c("SUBJID"))
        }
    }
} else {
    core_list <- c("c_keyvar","c_entry", "a_eendpt")

    # Too little pre-treatment entries:
    # "c_radio", "c_chemo","c_conmed"

    for (file in core_list) {
        imp <- read.csv(paste0(file_path, file, ".csv"))
        assign(strsplit(file, "_")[[1]][2], imp)
    }
}

keyvar <- dplyr::select(keyvar, all_of(c("SUBJID", "TXGROUP", "AGE", "SEXCD", "RACECD", "B_HGB", "B_SEREPO", "TUMORCD", "B_ECOGN", "B_LDHN")))
entry <- dplyr::select(entry, all_of(c("SUBJID", "ELIGYN")))
endpt <- dplyr::select(eendpt, all_of(c("SUBJID", "PFSCD", "PFSDY")))

X <- left_join(keyvar, entry, by = c("SUBJID"))
X <- left_join(X, endpt, by = c("SUBJID"))

# remove non-eligible patients
X <- filter(X, ELIGYN == 1)

pfs <- data.frame(T = X$PFSDY, C = X$PFSCD)
tx <- X$TXGROUP

X$PFSDY <- NULL
X$PFSCD <- NULL
X$TXGROUP <- NULL

X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

# No missing indicators to be corrected

Y_list <- impute_survival(T = pfs$T, C = pfs$C, X = X_imp)

NCT00119613 <- list(X_imp, Y_list, W)