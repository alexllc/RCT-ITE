
library(tidyverse)
library(broom)
library(ggfortify)

check_df <- data.frame(outcome = Y, trainX = X)
model <- lm(outcome ~ ., data = check_df)
vif_res <- car::vif(model)
model.diag.metrics <- augment(model)

pdf("NCT00364013_lm_plot.pdf")
par(mfrow = c(2, 2))
plot(model)
dev.off()

pdf("NCT00364013_diagnostic_plot.pdf")
autoplot(model)
dev.off()

model.diag.metrics <- model.diag.metrics %>% mutate(index = 1:nrow(model.diag.metrics)) 
  
# Inspect the data
head(model.diag.metrics)

# outliers
outlier_id <- filter(model.diag.metrics, abs(.std.resid) >= 3)
outlier_id <- outlier_id$index

# high leverage
high_lev_thresh <- 2*(dim(X)[2] + 1)/ dim(X)[1]
highlev_id <- filter(model.diag.metrics, .hat > high_lev_thresh)
highlev_id <- highlev_id$index


# high influence
high_influ_thresh <- 4/(dim(X)[1] - dim(X)[2] - 1)
highinflu_id <- filter(model.diag.metrics, .cooksd > high_influ_thresh)
highinflu_id <- highinflu_id$index

rmv_id <- unique(c(outlier_id, highlev_id, highinflu_id))

pdf("influence.pdf")
# Cook's distance
plot(model, 4, id.n = 5)
# Residuals vs Leverage
plot(model, 5)
dev.off()
