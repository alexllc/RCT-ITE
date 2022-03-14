library(hettx)

source("./audit/model_eval/load_lib.R")

clin_df <- read.csv("./dat/NCT00364013_adls.csv")

clin_df$PFS_mo <- clin_df$PFSDYCR / 30.4167
clin_df$OS_mo <- clin_df$DTHDY / 30.4167

X <- select(clin_df, all_of(c("LIVERMET", "DIAGMONS", "AGE",  "SEX", "B_WEIGHT", "B_HEIGHT", "RACE",  "B_ECOG", "HISSUBTY", "B_METANM", "DIAGTYPE")))
W <- clin_df$ATRT

########################################################################
#### Step 2: address outcome censoring
########################################################################

# 1. Efron's tail correction
source("./audit/original_surv.r")
efron_tail <- exp(impute.survival(surv.time = clin_df$PFS_mo, censor = clin_df$PFSCR)) # no need to convert to log scale

# change names from source data
X <- as.matrix(X)
Y <- efron_tail

best_tau <- tau_ptof_pred

test_hte <- data.frame(Y  = efron_tail, Z = W)

# enumerate our own all possible treatment assignments with the same proportion as the observed values
B <- 500
Z_tilde <- matrix(0, nrow = dim(clin_df)[1], ncol = B)
tau_tilde <- matrix(0, nrow = dim(clin_df)[1], ncol = B)

for (i in 1:B) {
    if (i == B/2) print("Halfway through.")

    Z_tilde[,i] <- sample(W)
    ptof <- PTOforest(x = X, tx = Z_tilde[,i], y = Y)
    tau_tilde[,i] <- predict(ptof, X)
}

hte_no_adj <- detect_idiosyncratic(formula = Y ~ Z, data = test_hte, plugin = TRUE, test.stat = "SKS.stat",  tau.hat = best_tau, te.vec = tau_tilde, n.cores = 8)

test_hte_adj <- cbind(test_hte, X)

ctrl_fml <- ~ DIAGMONS + AGE + SEX + B_WEIGHT + B_HEIGHT + RACE + B_ECOG
int_fml <- ~ LIVERMET + HISSUBTY + B_METANM + DIAGTYPE

hte_adj <- detect_idiosyncratic(formula = Y ~ Z, data = test_hte_adj, interaction.formula = int_fml, control.formula = ctrl_fml, test.stat = "SKS.stat.int.cov", plugin = FALSE, n.cores = 4)