library(rpart)
library(readr)
library(caTools)
library(dplyr)
library(party)
library(partykit)
library(rpart.plot)
library(C50)
library(printr)
library(RWeka)
library(FSelector)
library(randomForest)
library(reprtree)
library(e1071)
library(caret)

# Recovered vs Fatal outcome
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_vsurf.csv")
mydata <- multiplex[,-1]
mydata2 <- mydata[,1:85]
mydata2 <- as.data.frame(mydata2)

#select_cols <- colnames(multiplex)
#select_cols <- select_cols[!select_cols %in% c("Age", "BMI", "Weight","Temperature", "HR",	"RR",	"SpO2",	"Disease_Score", "Outcome", "Outcome2", "Category2", "Day")]
#multiplex <- multiplex[,select_cols]

mydata2$Outcome <- as.factor(mydata2$Outcome)

# Random Forests
# train model
set.seed(123)
sample_data = sample.split(mydata2, SplitRatio = 0.7)
train_data <- subset(mydata2, sample_data == TRUE)
test_data <- subset(mydata2, sample_data == FALSE)

set.seed(1)
rf_train <- randomForest(Outcome~., data=train_data, importance=T, ntree=1000, nodesize=5)
print(rf_train)
plot(rf_train)
reprtree:::plot.getTree(rf_train)
# summarize the fit
summary(rf_train)

# importance of each predictor
important <- importance(rf_train)
Important_Features <- data.frame(Feature = row.names(important), Importance = important [,3]) #to plot the Mean Decrease Accuracy
#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               Importance) , y = Importance) ) +
  geom_bar(stat = "identity", 
           fill = "#800080") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy")+
  ggtitle("Important Features in Random Forest\n") +
  theme(axis.text.y = element_text(size = 6))
plot_

#plot important variables
result <- data.frame(importance(rf_train,type=1))
result1 <- row.names(result)[result$MeanDecreaseAccuracy>3]
varImpPlot(rf_train,n.var=min(length(result1), nrow(rf_train$importance))) #same plot as above

varUsed(randomForest(Outcome~., data=train_data, tree=1000, count=T, by.tree=T))
plot(randomForest(Outcome~., data=train_data, importance=T, ntree=1000, nodesize=5), log="y")

# make predictions in the test data set
predictions <- predict(rf_train, test_data)
# summarize accuracy and predicitons in the test data set
table(predictions, test_data$Outcome)
mean(predictions == test_data$Outcome)
confusionMatrix(predictions, test_data$Outcome)

# other type of predictions in the test data set
predictions2 <- predict(rf_train, test_data, type="response", predict.all=TRUE)
table(predictions2, test_data$Outcome)
mean(predictions2 == test_data$Outcome)
confusionMatrix(predictions2, test_data$Outcome)

#performance of the random forests model applied
library(ROCR)
pred1=predict(rf_train,type = "prob")
pred = prediction(pred1[,2], train_data$Outcome)
# 1. True Positive and Negative Rate
pred2 = performance(pred, "tpr","fpr")
# 2. Plot the ROC curve
plot(pred2,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 3. Area under curve
auc <- performance(pred, measure = "auc")
auc
auc@y.values
# 4 - Accuracy
acc.perf = performance(pred, measure = "acc")
plot(acc.perf)

# to extract the maximum accuracy and the cutoff 
ind = which.max( slot(acc.perf, "y.values")[[1]] )
acc = slot(acc.perf, "y.values")[[1]][ind]
cutoff = slot(acc.perf, "x.values")[[1]][ind]
print(c(accuracy= acc, cutoff = cutoff))

#cut-off, sensitivity and specificity
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(pred2, pred))

# Train the RF model

set.seed(123)
sample_data = sample.split(mydata2, SplitRatio = 0.7) #split the dataset into train and test set in the ratio 70:30.
train_data <- subset(mydata2, sample_data == TRUE)
test_data <- subset(mydata2, sample_data == FALSE)
summary(train_data)
summary(test_data)

# create a Random Forest model with default parameters and then we will fine tune the model by changing ‘mtry’. 
# We can tune the random forest model by changing the number of trees (ntree) and 
# the number of variables randomly sampled at each stage (mtry).
rf_train2 <- randomForest(Outcome~., data=train_data, importance=T, ntree=1000, nodesize=5, maxnodes = NULL)
rf_train2

# Predicting on train set
predTrain <- predict(rf_train2, train_data, type = "class")

# Checking classification accuracy
table(predTrain, train_data$Outcome)  
mean(predTrain == train_data$Outcome)  

# Predicting on test set
predValid <- predict(rf_train2, test_data, type = "class")

# Checking classification accuracy
mean(predValid == test_data$Outcome)                    
table(predValid,test_data$Outcome)

# importance of each predictor
important <- importance(rf_train2)
Important_Features <- data.frame(Feature = row.names(important), Importance = important [,3])
#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               Importance) , y = Importance) ) +
  geom_bar(stat = "identity", 
           fill = "#800080") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy")+
  ggtitle("Important Features in Random Forest\n") +
  theme(axis.text.y = element_text(size = 6))
plot_ #same results as the model in tutorial 1

#performance of the random forests model applied
library(ROCR)
pred1=predict(rf_train2,type = "prob")
pred = prediction(pred1[,2], train_data$Outcome)
# 1. True Positive and Negative Rate
pred2 = performance(pred, "tpr","fpr")
# 2. Plot the ROC curve
plot(pred2,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 3. Area under curve
auc <- performance(pred, measure = "auc")
auc
auc@y.values
# 4 - Accuracy
acc.perf = performance(pred, measure = "acc")
plot(acc.perf)

# to extract the maximum accuracy and the cutoff 
ind = which.max( slot(acc.perf, "y.values")[[1]] )
acc = slot(acc.perf, "y.values")[[1]][ind]
cutoff = slot(acc.perf, "x.values")[[1]][ind]
print(c(accuracy= acc, cutoff = cutoff))

#cut-off, sensitivity and specificity
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(pred2, pred))

saveRDS(rf_train,"rftrain_model1.rds")
saveRDS(rf_train2,"rftrain_model2.rds")

rf_train <- readRDS("rftrain_model1.rds")
rf_train2 <- readRDS("rftrain_model2.rds")

# Tuning and training RF model

# Tuning a model is very tedious work. There are lot of combination possible between the parameters. 
# You don't necessarily have the time to try all of them. 
# A good alternative is to let the machine find the best combination for you. There are two methods available:
# Grid search and Random search

# Set the control parameter: You will proceed as follow to construct and evaluate the model:
# Evaluate the model with the default setting
# Find the best number of mtry
# Find the best number of maxnodes
# Find the best number of ntrees
# Evaluate the model on the test dataset

# Step 1- K-fold cross validation is controlled by the trainControl() function
# Define the control
trControl <- trainControl(method = "cv",
                          number = 10,   #10-fold cross validation
                          search = "grid")
# You will use caret library to evaluate your model. The library has one function called train() to evaluate machine learning algorithms.
rf_default <- train(Outcome~., data=train_data, method = "rf", metric= "Accuracy", trControl = trainControl(), tuneGrid = NULL)
print(rf_default)
# choose mtry with the higher accuracy score

#Step 2- seach best mtry
set.seed(1234)
tuneGrid <- expand.grid(.mtry = c(1: 10))
rf_mtry <- train(Outcome~.,
                 data = train_data,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 5,
                 ntree = 300)
print(rf_mtry)
# Accuracy was used to select the optimal model using  the largest value.
#The best value of mtry is stored in:
rf_mtry$bestTune$mtry
#You can store it and use it when you need to tune the other parameters.
max(rf_mtry$results$Accuracy)
best_mtry <- rf_mtry$bestTune$mtry 
best_mtry

#Step 3 - Search the best maxnodes
#You need to create a loop to evaluate the different values of maxnodes. In the following code, you will:
#Create a list
#Create a variable with the best value of the parameter mtry; Compulsory
#Create the loop
#Store the current value of maxnode
#Summarize the results
store_maxnode <- list()  #The results of the model will be stored in this list
tuneGrid <- expand.grid(.mtry = best_mtry) #Use the best value of mtry
for (maxnodes in c(15: 30)) {   #Compute the model with values of maxnodes starting from 5 to 15.
  set.seed(1234)
  rf_maxnode <- train(Outcome~.,
                      data = train_data,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 5,
                      maxnodes = maxnodes,   #For each iteration, maxnodes is equal to the current value of maxnodes
                      ntree = 300)
  current_iteration <- toString(maxnodes)    # Store as a string variable the value of maxnode.
  store_maxnode[[current_iteration]] <- rf_maxnode    #Save the result of the model in the list.
}
results_mtry <- resamples(store_maxnode)  #Arrange the results of the model
summary(results_mtry)  #Print the summary of all the combination.
#The last value of maxnode has the highest accuracy. You can try with higher values to see if you can get a higher score.

#Step 4 - Search the best ntrees
#Now that you have the best value of mtry and maxnode, you can tune the number of trees. The method is exactly the same as maxnode.
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(5678)
  rf_maxtrees <- train(Outcome~.,
                       data = train_data,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 5,
                       maxnodes = 22,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
#Record the OOB error rate and see the number of trees where the out of bag error rate stabilizes and reach minimum.

fit_rf <- train(Outcome~.,
                train_data,
                method = "rf",
                metric = "Accuracy",
                tuneGrid = tuneGrid, #best mtry value is here
                trControl = trControl,
                importance = TRUE,
                nodesize = 5, #change this if necessary
                ntree = 2000,   #change this if necessary
                maxnodes = 15)  #change this if necessary

# Step 5 - Evaluate the model in the test data
# The library caret has a function to make prediction.
prediction <-predict(fit_rf, test_data)
table(prediction, test_data$Outcome)
mean(prediction == test_data$Outcome)
# You can use the prediction to compute the confusion matrix and see the accuracy score
confusionMatrix(prediction, test_data$Outcome)


# Step 6 - Visualize Result
varImpPlot(fit_rf) #This function only works for objects of class `randomForest'
#Important features are likely to appear closer to the root of the tree, while less important features will often appear closed to the leaves.

#Step 7 - running the rf model after setting up the parameters 
set.seed(1)
rf_train3 <- randomForest(Outcome~., data=train_data, importance=T, mtry=2, maxnodes = 15, ntree=1000, nodesize=5)
rf_train3
plot(rf_train3)
reprtree:::plot.getTree(rf_train3)

# summarize the fit
summary(rf_train3)

saveRDS(rf_train3,"rftrain_model3.rds")
rf_train3 <- readRDS("rftrain_model3.rds")

# Predicting on train set
predTrain <- predict(rf_train3, train_data, type = "class")
# Checking classification accuracy
table(predTrain, train_data$Outcome)  
mean(predTrain == train_data$Outcome)  
confusionMatrix(predTrain, train_data$Outcome)

# Predicting on test set
predValid <- predict(rf_train3, test_data, type = "class")
# Checking classification accuracy
mean(predValid == test_data$Outcome)                    
table(predValid,test_data$Outcome)
confusionMatrix(predValid, test_data$Outcome)


# importance of each predictor
important <- importance(rf_train3, type = 1) #type = 1, to extract only the Mean Decreased Accuracy results
Important_Features <- data.frame(Feature = row.names(important), Importance = important [,1])
#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               Importance) , y = Importance) ) +
  geom_bar(stat = "identity", 
           fill = "#800080") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy")+
  ggtitle("Important Features in Random Forest\n") +
  theme(axis.text.y = element_text(size = 6))
plot_
  
#plot important variables
result <- data.frame(importance(rf_train3,type=1))
result1 <- row.names(result)[result$MeanDecreaseAccuracy>3]
varImpPlot(rf_train3,n.var=min(length(result1), nrow(rf_train3$importance)))

#performance of the random forests model applied
library(ROCR)
pred1=predict(rf_train3,type = "prob")
pred = prediction(pred1[,2], train_data$Outcome)
# 1. True Positive and Negative Rate
pred2 = performance(pred, "tpr","fpr")
# 2. Plot the ROC curve
plot(pred2,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
# 3. Area under curve
auc <- performance(pred, measure = "auc")
auc
auc@y.values
# 4-  Accuracy
acc.perf = performance(pred, measure = "acc")
plot(acc.perf)
 
# to extract the maximum accuracy and the cutoff 
ind = which.max( slot(acc.perf, "y.values")[[1]] )
acc = slot(acc.perf, "y.values")[[1]][ind]
cutoff = slot(acc.perf, "x.values")[[1]][ind]
print(c(accuracy= acc, cutoff = cutoff))

#cut-off, sensitivity and specificity
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(pred2, pred))

# Model 4 - Using C5.0 trees
# fit model
fit <- C5.0(Outcome~., data=mydata2, trials=100)
# summarize the fit
print(fit)
# make predictions
predictions <- predict(fit, mydata2)
# summarize accuracy
table(predictions, mydata2$Outcome)
plot(fit)

#Training and testing the model
set.seed(123)
sample_data = sample.split(mydata2, SplitRatio = 0.7)
train_data <- subset(mydata2, sample_data == TRUE)
test_data <- subset(mydata2, sample_data == FALSE)

model <- C5.0(Outcome~., data=train_data, trials=100)
print(model)
results <- predict(object=model, newdata=test_data, type="class")
table(results, test_data$Outcome)
plot(model)


# Model 5 - Using C4.5 trees
m1 <- J48(Outcome~., data = mydata2)
if(require("party", quietly = TRUE)) plot(m1)
summary(m1)
# make predictions
predictions <- predict(m1, mydata2)
# summarize accuracy
table(predictions, mydata2$Outcome)

library(FSelector)
information.gain(Outcome~., data = mydata2)

#First, use caret to create a 10-fold training set. Then train the model.
library(RWeka)
library(caret)
library(e1071)

set.seed(123)
sample_data = sample.split(mydata2, SplitRatio = 0.7)
train_data <- subset(mydata2, sample_data == TRUE)
test_data <- subset(mydata2, sample_data == FALSE)

set.seed(1958)  # set a seed to get replicable results
train <- createFolds(train_data$Outcome, k=10)
C45Fit <- train(Outcome~., method="J48", data=train_data,
                tuneLength = 5,
                trControl = trainControl(
                  method="cv", indexOut=train))
#cv is the K-fold cross-validation algorithm
#K-fold cross-validation (CV) is a robust method for estimating the accuracy of a model.
#A less obvious but potentially more important advantage of k-fold CV is that it often gives more accurate estimates of the test error rate than does LOOCV (James et al. 2014).
#In practice, one typically performs k-fold cross-validation using k = 5 or k = 10, 
#as these values have been shown empirically to yield test error rate estimates that suffer neither from excessively high bias nor from very high variance.

summary(C45Fit)
C45Fit
C45Fit$finalModel
plot(C45Fit)
plot(C45Fit$finalModel)

# Model 6 - Using e1071
library(e1071) #Author DataFlair
svm_fit <- svm(Outcome~., data = mydata2)
plot(svm_fit, data = mydata2)
pred <- predict(model_svm,x)
confusionMatrix(pred,y$Species)

# Model 7 - Using CART trees
library(rpart)
library(rpart.plot)
library(rattle)
# Create a decision tree model
set.seed(123)
sample_data = sample.split(mydata2, SplitRatio = 0.7)
train_data <- subset(mydata2, sample_data == TRUE)
test_data <- subset(mydata2, sample_data == FALSE)

model <- rpart(Outcome~., data=train_data, method = "class", minsplit=15, cp=-1)
print(model)

rpart.plot(model, box.palette="RdBu", shadow.col="gray", nn=TRUE)
fancyRpartPlot(model)
plotcp(model)
printcp(model)
rsq.rpart(model)
summary(model)

# Model 8 - ctree
#conditional parting plot as follows:
model2 <- ctree(Outcome~., data=train_data)
plot(model2)

# make predictions
results <- predict(object=model, newdata=test_data, type="class")
# summarize accuracy
table(results, test_data$Outcome)

#prune the tree
#Prune back the tree to avoid overfitting the data. 
#Typically, you will want to select a tree size that minimizes the cross-validated error, 
#the xerror column printed by printcp( ).
#Specifically, use printcp( ) to examine the cross-validated error results, 
#select the complexity parameter associated with minimum error, 
#and place it into the prune( ) function. Alternatively, you can use the code fragment

fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]
#to automatically select the complexity parameter associated with the smallest cross-validated error.

tree2 <- prune(model, cp=-1)
# plot the pruned tree 
rpart.plot(tree2, box.palette="RdBu", shadow.col="gray", nn=TRUE)
plotcp(tree2)
printcp(tree2)
rsq.rpart(tree2)
summary(tree2)

tree <- rpart(Outcome~., data=test_data, method = "class", minsplit=25, cp=-1)
# minsplit = 0, minbucket = 0
# Visualize the decision tree with rpart.plot
rpart.plot(tree, box.palette="RdBu", shadow.col="gray", nn=TRUE)
fancyRpartPlot(tree)
plotcp(tree)
printcp(tree)
rsq.rpart(tree)
summary(tree)

#prune the tree
tree2 <- prune(tree, cp=-1)
rpart.plot(tree2, box.palette="RdBu", shadow.col="gray", nn=TRUE)
lotcp(tree2)
printcp(tree2)
rsq.rpart(tree2)
summary(tree2)




