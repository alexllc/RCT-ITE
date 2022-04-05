source("./bin/load_lib.R")
source("./bin/impute_survival.R")

file_path <- "./dat/PDS/LungSm_Amgen_2002_266_NCT00119613/csv/"

keyvar <- read.csv(paste0(file_path, "c_keyvar.csv"))
bchar <- read.csv(paste0(file_path, "c_bchar.csv"))
medhis <- read.csv(paste0(file_path, "c_medhis.csv"))
radio <- read.csv(paste0(file_path, "c_radio.csv"))
c_conmed
c_vitals
c_hemat # (pretreatment)
c_chem

# No missing indicators to be corrected

keyvar_sel <- dplyr::select(eval, all_of(c("SUBJID", "AGE", "SEXCD", "RACECD", "B_HGB", "B_SEREPO", "TUMORCD", "B_ECOGN", "B_LDHN")))

bchar_sel <- dplyr::select(eval, all_of(c("SUBJID", "B_WEIGHT", "B_HEIGHT")))

medhis_sel <- dplyr::select(eval, all_of(c("SUBJID", "MEDHXCD", "MHSTATCD")))

radio_sel <- dplyr::select(eval, all_of(c("SUBJID", )))

X <- dplyr::select(eval, all_of(c()))