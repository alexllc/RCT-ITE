source("/home/alex/Documents/lab/RCT-ITE/bin/impute_survival.R")


corevar <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_263_NCT00339183/csv/corevar.csv")
eendpt <- read.csv("/home/alex/Documents/lab/RCT-ITE/dat/PDS/Colorec_Amgen_2006_263_NCT00339183/csv/a_eendpt.csv")

corevar[corevar == "" | corevar == "NA"] <- NA

X <- dplyr::select(corevar, all_of(c("AGE", "SEXCD" , "RACCATCD", "PRBEVCDR", "PRBEVR"  , "PRBEVCDM", "PRBEVM" , "PRBEVCD" , "PRBEV", "PROXACDR", "PROXAR"  , "PROXACDM", "PROXAM"  , "PROXALCD", "PROXAL"  , "B_METANM", "LIVRONCD", "PADJYNCD" , "DPADJUV" , "B_LDHYN" , "B_LDH2YN", "B_LDHNM" , "BECOGICD", "B_ECOGI" , "DIAGTYPE", "TUMCAT"  , "LIVERMET", "KRASCD", "B_ECOGCD", "TRTDUR")))
X <- missing_too_much(X)


W <- as.numeric(corevar$ATRT == "Panitumumab + FOLFIRI")

X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

Y_list <- impute_survival(T = eendpt$PFSDYCR / 30.4167, C = eendpt$PFSCR, X = X_imp)

NCT00339183 <- list(X_imp, Y_list, W)