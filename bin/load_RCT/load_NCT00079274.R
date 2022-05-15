#' Script to load data from the trial NCT00079274
#' Arms B, C, E, and F were discontinued as of June 1, 2005.

file_path <- "./dat/PDS/Colorec_Allianc_2004_161_NCT00079274/"

char <- read.csv(paste0(file_path, "characteristic.csv"))
objectives <- read.csv(paste0(file_path, "objectives.csv"))
tox <- read.csv(paste0(file_path, "tox.csv"))

rand_char <- filter(char, (ARM == "A" | ARM == "D") & is.na(EXCLUDED))
rand_obj <- objectives[objectives$mask_id %in% rand_char$mask_id, ]
rand_tox <- tox[tox$mask_id %in% rand_char$mask_id, ]

X <- dplyr::select(rand_char, all_of(c("agecat", "racecat", "SEX", "BWL_OBS", "BWL_PERF", "HISTO_G", "NODES", "NUMCYCLE", "STAGE_G", "PS", "wild", "bmi2")))
X <- missing_too_much(X)
X_imp <- impute_df_missing(clin_df = X, save_ddt = FALSE)

W <- as.numeric(rand_char$ARM == "A")

NCT00079274<- list(X_imp, W)

# In May 2004, the Oncologic Drugs Advisory Committee to the FDA
# recommended that disease-free survival become an acceptable regulatory
# endpoint. This was based on data presented from a large pooled analysis
# showing a very tight relationship between disease-free and overall
# survival in previous trials in adjuvant colon cancer. Accordingly, we are
# changing the primary endpoint of our trial to disease-free survival, which
# will allow the primary analysis to occur one year sooner than with the
# overall survival endpoint.

toxicity <- rand_tox %>% group_by(mask_id) %>% mutate(total_ae_grade = sum(GRADE)) %>% select(all_of(c("mask_id", "total_ae_grade"))) %>% unique()
toxicity <- toxicity[toxicity$mask_id %in% rand_char$mask_id,]

outcomes <- left_join(rand_obj, toxicity, by = c("mask_id"))
outcomes$total_ae_grade[is.na(outcomes$total_ae_grade)] <- 0

OS <- data.frame(T = outcomes$futime8 / 30.4167, C = outcomes$fustat8)
DFS <- data.frame(T = outcomes$dfstime5 / 30.4167, C = outcomes$dfsstat5)
PFS <- data.frame(T = outcomes$pgtime5 / 30.4167, C = outcomes$pgstat5)
AE <- outcomes$total_ae_grade

NCT00079274_outcomes <- c("OS", "DFS", "PFS", "AE")
for (outcome in NCT00079274_outcomes) {
    if (outcome != "AE") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp)))
    } else {
        assign(paste0(outcome, "_Y_list"), list(outcomes[["total_ae_grade"]], NA, NA))
    }
}

