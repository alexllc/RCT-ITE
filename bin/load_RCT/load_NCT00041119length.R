source("./bin/load_RCT/load_NCT00041119chemo.R")

W <- as.numeric(eval$indrx == 2 | eval$indrx == 4) # Experimental arm is: 2=CA-6 or 4=T-6
chemo_type <- NCT00041119chemo[[2]] # CA = 1 vs Paclitaxel = 0
X_imp <- cbind(X_imp, chemo_type)

NCT00041119length <- list(X_imp, W)

NCT00041119length_outcomes <- c("DFS", "OS", "AE")

for (outcome in NCT00041119length_outcomes) {
    if (outcome != "AE") {
        assign(paste0(outcome, "_Y_list"), do.call(impute_survival, list(T = get(outcome)[,1], C = ceiling(get(outcome)[,2]), X = X_imp)))
    } else {
        assign(paste0(outcome, "_Y_list"), list(get(outcome), NA, NA))
    }
}

save(NCT00041119length, NCT00041119length_outcomes, OS_Y_list, DFS_Y_list, AE_Y_list, file = "./bin/load_RCT/RCT_obj/NCT00041119length.RData")