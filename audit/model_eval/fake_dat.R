library(ganGenerativeData)

source("./audit/model_eval/load_data.R")

original_df <- cbind(efron_tail_Yn, tx, as.matrix(covar))
gan(original_df, epochs = 2000)
