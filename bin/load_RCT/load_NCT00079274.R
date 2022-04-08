source("./bin/impute_survival.R")

clinchar <- read.csv("./dat/PDS/Colorec_Allianc_2004_161_NCT00079274/characteristic.csv")

objectives <- read.csv("./dat/PDS/Colorec_Allianc_2004_161_NCT00079274/objectives.csv")

rand_char <- filter(clinchar, (ARM == "A" | ARM == "D") & is.na(EXCLUDED))
rand_obj <- objectives[objectives$mask_id %in% rand_char$mask_id, ]

X <- dplyr::select(rand_char, all_of(c("agecat", "racecat", "SEX", "BWL_OBS", "BWL_PERF", "HISTO_G", "NODES", "STAGE_G", "PS", "wild", "bmi2")))

x <- missing_too_much(X)

for(col in colnames(X)) {
    
}

W <- as.numeric(rand_char$ARM == "A")
table(W)

X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

# In May 2004, the Oncologic Drugs Advisory Committee to the FDA
# recommended that disease-free survival become an acceptable regulatory
# endpoint. This was based on data presented from a large pooled analysis
# showing a very tight relationship between disease-free and overall
# survival in previous trials in adjuvant colon cancer. Accordingly, we are
# changing the primary endpoint of our trial to disease-free survival, which
# will allow the primary analysis to occur one year sooner than with the
# overall survival endpoint.

Y_list <- impute_survival(T = rand_obj$dfstime5 / 30.4167, C = rand_obj$dfsstat5, X = X_imp)

NCT00079274<- list(X_imp, Y_list, W)