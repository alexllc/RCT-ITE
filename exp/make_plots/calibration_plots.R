source("./bin/load_lib.R")
library("RColorBrewer")

#
# Summarize existing calib
#
reest_trials <- c("NCT00003299", "NCT00041119chemo", "NCT00113763", "NCT00115765", "NCT00119613")

calib_raw_files <- list.files("./res/calib/backup/")
calib_raw_trials <- unique(unlist(lapply(strsplit(calib_raw_files, "_|\\."), function(X) X[1])))
calib_raw_trials <- calib_raw_trials[-1]

calib_raw_trials <- calib_raw_trials[!(calib_raw_trials %in% reest_trials)]


for (trial in calib_raw_trials) {
    if (!file.exists(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))) {
        source(paste0("./bin/load_RCT/load_", trial, ".R"))
    } else {
        load(paste0("./bin/load_RCT/RCT_obj/", trial, ".RData"))
    }

    trial_calib_list <- list()

    trial_files <- calib_raw_files[grep(paste0("^", trial,"*"), calib_raw_files)]
    outcome_list <- unique(unlist(lapply(strsplit(trial_files, "_|\\."), function(X) X[5])))

    for (outcome in outcome_list) {
        Y_list <- get(paste0(outcome, "_Y_list"))

        imp_type_Y <- "efronYn"
        Y <- as.numeric(Y_list[[2]]) # use Efron+Yn if available

        if (any(is.na(Y))) { # either does not have largest censroed outcome or it is not a time to event type outcome
            Y <- as.numeric(Y_list[[1]])
            imp_type_Y <- "efron"
        }

        outcome_files <- trial_files[grep(paste0("_", outcome, "_"), trial_files)]
        tau_method_used <- unique(unlist(lapply(strsplit(outcome_files, "_|\\."), function(X) X[2])))

        for (tau_meth in tau_method_used) {
            iter_files <- trial_files[grep(paste0(tau_meth, "_calib_", imp_type_Y, "_", outcome,"_"), trial_files)]
                iter_list <- unique(unlist(lapply(strsplit(iter_files, "_|\\."), function(X) X[6])))
                avg_calib_list <- list()
                for (iter in iter_list) {
                    calib <- read.csv(paste0("./res/calib/backup/", trial, "_", tau_meth, "_calib_", imp_type_Y, "_", outcome, "_", iter, "_res.csv"))
                    avg_calib_list <- append(avg_calib_list, list(c(calib[1,], calib[2,])))
                }
                avg_calib_sum <- do.call(rbind.data.frame, avg_calib_list)
                avg_calib_sum$term <- NULL
                avg_calib_sum$`term.1` <- NULL

                trial_calib_list <- append(trial_calib_list, list(c(outcome, tau_meth, as.list(colMeans(avg_calib_sum, na.rm = TRUE)))))
        }

    }
    trial_avg_calib <- do.call(rbind.data.frame, trial_calib_list)
    colnames(trial_avg_calib) <- c("outcome", "tau_method", "ate.estimate", "ate.std.error", "ate.statistic", "ate.pval", "hte.estimate", "hte.std.error", "hte.statistic", "hte.pval")

    write.csv(trial_avg_calib, paste0("./res/calib/backup/calib_avg/", trial, "_calib_res.csv"), row.names = FALSE)
}





#
# make calibration plots with summarized tables
#

calib_file_path <- list.files("./res/calib/backup/calib_avg/")
calib_trials <- unique(unlist(lapply(strsplit(calib_file_path, "_|\\."), function(X) X[1])))

for (trial in calib_trials) {
    calib_sum <- read.csv(paste0("./res/calib/backup/calib_avg/", trial, "_calib_res.csv"))

    calib_sum$p.sig <- rep("", dim(calib_sum)[1])
    calib_sum$p.sig[which(calib_sum$hte.pval <= 0.1)] <- "."
    calib_sum$p.sig[which(calib_sum$hte.pval <= 0.05)] <- "*"
    calib_sum$p.sig[which(calib_sum$hte.pval <= 0.01)] <- "**"
    calib_sum$p.sig[which(calib_sum$hte.pval <= 0.001)] <- "***"
    

    pd <- position_dodge(.7)    ### How much to jitter the points on the plot

    pdf(paste0("./res/calib/backup/calib_avg/", trial, "_calib_sum_plot.pdf"))
    print(ggplot(calib_sum,                ### The data frame to use.
    aes(x     = tau_method,
        y     = hte.estimate,
        color = outcome)) +

    geom_point(shape = 15,
            size  = 4,
            position = pd) +
    scale_color_brewer(palette = "Paired") +
    ylim(-2, 3) +

    geom_errorbar(aes(ymin  = hte.estimate - hte.std.error,
                    ymax  = hte.estimate + hte.std.error),
                    width = 0.1,
                    size  = 0.7,
                    position = pd) +
    geom_text(aes(x = tau_method, y = 2.8, label = p.sig), position = pd, size = 7) +
    geom_hline(yintercept=1, linetype="dashed", color = "#1F78B4") +
    theme(axis.title = element_text(face = "bold")) +
    ggtitle(trial) +
    ylab("Best linear fit ITE estimate") + 
    xlab("ITE estimator") + 
    theme_bw() +
    theme( plot.title = element_text(color="#6A3D9A", size=36, face="bold.italic", hjust = 0.5),
    axis.title.x = element_text(color="#33A02C", size=20, face="bold"),
    axis.title.y = element_text(color="#B2DF8A", size=20, face="bold"), 
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    )
    dev.off()
}
