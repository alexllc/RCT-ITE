library(tidyverse)
library(caret)

lm_df <- as.data.frame(cbind(OS_time = Y, X))
model <- lm(OS_time~., data = lm_df)
predictions <- model %>% predict(lm_df)
data.frame(RMSE = RMSE(predictions, lm_df$OS_time), R2 = caret::R2(predictions, lm_df$OS_time))

vif_res <- car::vif(model)
X_vif <- X[,vif_res < 2]

X <- X_vif