source("./bin/load_lib.R")

library(ComplexHeatmap)


biomarker_trials <- c("NCT00113763", "NCT00364013")

for (trial in biomarker_trials) {

    message(paste0(rep("=", 80)))
    message(paste0("Running trial: ", trial))
    message(paste0(rep("=", 80)))

    source(paste0("./bin/load_RCT/load_", trial, ".R"))
    X <- as.matrix(get(trial)[[1]])
    W <- get(trial)[[2]]
    outcome_list <- get(paste0(trial, "_outcomes"))

    tau_df <- data.frame(seq(1:dim(X)[1]))
    for (outcome in outcome_list) {
        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }
        tau <- try(read.csv(paste0("./res/ite_tau_estimates/backup/", trial, "_", imp_type_Y, "_", outcome, "_Rstack_tau_estimates.csv")))
        
        tau_df <- cbind(tau_df, tau[,1])

    }
    tau_df[,1] <- NULL

    scaled_tau <- as.data.frame(normalize.quantiles.robust(x = as.matrix(tau_df)))
    colnames(scaled_tau) <- outcome_list
    biomark <- as.data.frame(X[,grep("KRAS|NRAS|BRAF", colnames(X))])
    colnames(biomark) <- gsub("\\(|\\)", "", gsub("/", "_", colnames(biomark)))
    biomark <- as.data.frame(biomark)
    biomark[biomark == 1] <- "mutant"
    biomark[biomark==2] <- "wild-type"
    biomark[biomark != "mutant" & biomark != "wild-type"] <- "missing"
        

    col <- list(KRAS_exon_2_c12_13 = c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                KRAS_exon_3_c61 = c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                KRAS_exon_4_c117_146= c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                NRAS_exon_2_c12_13 = c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                NRAS_exon_4_c117_146 = c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                NRAS_exon_3_c61= c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                BRAF_exon_15_c600= c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                KRAS_exon_3_c59_c61= c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff"),
                NRAS_exon_3_c59_c61= c("mutant" = "#ff5900", "wild-type" = "#41d03f", "missing" = "#ffffff")
                )

    row_ha <- rowAnnotation(df = biomark, col = col)
    pdf(paste0(trial, "_heatmap.pdf"))
    Heatmap(scaled_tau, name = trial, right_annotation = row_ha, cluster_columns = FALSE)
    dev.off()

}