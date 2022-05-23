library(data.table)
library(rpart)
library(rpart.plot)
library(caret)
library(xgboost)
library(pROC)

set.seed(123)
full = fread('./HR_comma_sep.csv', stringsAsFactors = T)
full = full[sample(.N)]

#### Add Random Noise

tmp_std = sd(full[,satisfaction_level])
full[,satisfaction_level:=satisfaction_level + runif(.N,-tmp_std,tmp_std)]
full[,satisfaction_level:=satisfaction_level - min(satisfaction_level)]
full[,satisfaction_level:=satisfaction_level / max(satisfaction_level)]

tmp_std = sd(full[,last_evaluation])
full[,last_evaluation:=last_evaluation + runif(.N,-tmp_std,tmp_std) ]
full[,last_evaluation:=last_evaluation - min(last_evaluation)]
full[,last_evaluation:=last_evaluation / max(last_evaluation)]

tmp_min = min(full[,number_project])
tmp_std = sd(full[,number_project])
full[,number_project:=number_project + sample(-ceiling(tmp_std):ceiling(tmp_std),.N, replace=T)]
full[,number_project:=number_project - min(number_project) + tmp_min]

tmp_min = min(full[,average_montly_hours])
tmp_std = sd(full[,average_montly_hours])
full[,average_montly_hours:=average_montly_hours + sample(-ceiling(tmp_std):ceiling(tmp_std),.N, replace=T)]
full[,average_montly_hours:=average_montly_hours - min(average_montly_hours) + tmp_min]

tmp_min = min(full[,time_spend_company])
tmp_std = sd(full[,time_spend_company])
full[,time_spend_company:=time_spend_company + sample(-ceiling(tmp_std):ceiling(tmp_std),.N, replace=T)]
full[,time_spend_company:=time_spend_company - min(time_spend_company) + tmp_min]

tmp_min = min(full[,number_project])
tmp_std = sd(full[,number_project])
full[,number_project:=number_project + sample(-ceiling(tmp_std):ceiling(tmp_std),.N, replace=T)]
full[,number_project:=number_project - min(number_project) + tmp_min]

#### Create Train / Test and Folds

train = full[1:12000]
test = full[12001:14999]

cv <- createFolds(train[,left], k = 10)
# Control
ctrl <- trainControl(method = "cv",index = cv)#### Train Tree
tree.cv <- train(x = train[,-"left"], y = as.factor(train[,left]), method = "rpart2", tuneLength = 7,     trControl = ctrl, control = rpart.control())

tree.model = tree.cv$finalModel

rpart.plot(tree.model, type = 2,extra =  7,fallen.leaves = T)
rpart.plot(tree.model, type = 2,extra =  2,fallen.leaves = T)

tree.preds = predict(tree.model, test)[,2]
tree.roc_obj <- roc(test[,left], tree.preds)

cat("Tree AUC ", auc(tree.roc_obj))