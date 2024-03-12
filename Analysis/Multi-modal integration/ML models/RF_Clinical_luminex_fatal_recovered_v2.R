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
library(tidyverse)
library(tidymodels)

# Recovered vs Fatal outcome
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_vsurf.csv")
mydata <- multiplex[,-1]
mydata2 <- mydata[,1:85]
mydata2 <- as.data.frame(mydata2)

# retirando as colunas indesejadas para o heatmap
#select_cols <- colnames(multiplex)
#select_cols <- select_cols[!select_cols %in% c("Age", "BMI", "Weight","Temperature", "HR",	"RR",	"SpO2",	"Disease_Score", "Outcome", "Outcome2", "Category2", "Day")]
#multiplex <- multiplex[,select_cols]

mydata2$Outcome <- as.factor(mydata2$Outcome)

# Reference tutorial: https://notast.netlify.app/post/explaining-predictions-random-forest-post-hoc-analysis-permutation-impurity-variable-importance/
# Explaining Predictions: Random Forest Post-hoc Analysis (permutation & impurity variable importance)

# train/test set
set.seed(123)
sample_data = sample.split(mydata2, SplitRatio = 0.7)
train_data <- subset(mydata2, sample_data == TRUE)
test_data <- subset(mydata2, sample_data == FALSE)

# create recipe object
covid_recipe<-recipe(Outcome~., data=train_data) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal(), -Outcome)# need dummy variables to be created for some `randomForestexplainer` functions though random forest model itself doesnt need explicit one hot encoding 
# process the traing set/ prepare recipe(non-cv)
covid_prep <-covid_recipe %>% prep(training = train_data, retain = TRUE)

set.seed(123)
data_split <- initial_split(mydata2, prop=0.7, strata = "Outcome")
data_train <- training(data_split)
data_test <- testing(data_split)

# create recipe object 2
covid_recipe2<-recipe(Outcome~., data=data_train) %>%
  step_impute_knn(all_predictors()) %>% 
  step_dummy(all_nominal(), -Outcome)# need dummy variables to be created for some `randomForestexplainer` functions though random forest model itself doesnt need explicit one hot encoding 
# process the traing set/ prepare recipe(non-cv)
covid_prep2 <-covid_recipe2 %>% prep(training = data_train, retain = TRUE)

# 1- Random forest

# 1.1- Ensemble Models (tree-based)

# Random forests and gradient boosting are also examples of tree based ensemble models. 
# In tree based ensemble models, multiple variants of a simple tree are combined to develop a more complex model with higher prediction power. 
# The trade off of the ensemble is that the model becomes a black box and cannot be interpreted as simply as a decision tree. 
# Thus, post-hoc analysis is required to explain the ensemble model.

# 1.2- Random forest mechanism

# Random forest is a step up of bootstrap aggregating/ bagging and bootstrap aggregating is a step up of decision tree. 
# In bagging, multiple bootstrap samples are of the training data are used to train multiple single regression tree. 
# As bootstrap samples are used, the same observation from the training set can be used more than once to train a regression tree. 
# For regression problems, the ensemble occurs when the predictions of each tree are averaged to provide the predicted value for the bagging model. 
# For classification problems, the ensemble occurs when the class probabilities of each tree are averaged to provide the predicted class. 
# Alternatively, the majority predicted class of the ensemble of trees will be the predicted class for the model. 
# Random forest is similar to bootstrap aggregating except that a random subset of variables is used for each split. 
# The number of variables in this random subset of variables is known as mtry. The default mtry in regression problems is a third of the total number of available variables in the dataset. 
# The default mtry in classification problems is the square root of the total number of available variables. - that is why the RF models 1, 2 and 3 and VSURF use 9 variables/tree.

# 1.3- Training random forest models
# We will train two random forest where each model adopts a different ranking approach for feature importance. The two ranking measurements are:

# 1.3.1- Permutation based.
# Permuting values in a variable decouples any relationship between the predictor and the outcome which renders the variable pseudo present in the model. 
# When the out of bag sample is passed down the model with the permutated variable, the increased in error/ decreased in accuracy of the model can be calculated. 
# This accuracy/error is compared against the baseline accuracy/error (for regression problems, the measurement will be R2) to calculate the importance score. 
# This computation is repeated for all variables. 
# The variable with the largest decrease in accuracy/largest increase in error is considered the most important variable because not considering the variable has the largest penalty on prediction. 
# In the case of the randomForest package, the score measures the decrease in accuracy and this score is refer as the MeanDecreaseAccuracy.

# Impurity based. Impurity refers to gini impurity/ gini index. 
# The concept of impurity for random forest is the same as regression tree. 
# Features which are more important have a lower impurity score/ higher purity score/ higher decrease in impurity score. 
# The randomForest package, adopts the latter score which known as MeanDecreaseGini.
  
set.seed(1)
rf_model<-rand_forest(trees = 1000,  mode = "classification") %>% 
  set_engine("randomForest",  importance=T, localImp = T, ) %>% # importance = T to have permutation score calculated
  fit(Outcome~ ., data = juice(covid_prep)) # localImp=T for randomForestExplainer(next post)

saveRDS(rf_model,"rftrain_model4.rds")
rf_model <- readRDS("rftrain_model4.rds")

# 2- Feature Importance

# We can extract permutation based importance as well as gini based importance directly from the model trained by the randomForest package.

# 2.1- Permutation-based importance

# Using the tidyverse approach to the extract results, remember to convert MeanDecreaseAccuracy from character to numeric form for arrange to sort the variables correctly. 
# Otherwise, R will recognise the value based on the first digit while ignoring log/exp values. 
# For instance, if MeanDecreaseAccuracy was in character format, rest_ecg_ST.T.abnorm will have the highest value as R recognises the first digit in 4.1293934444743e-05 and fails to recognise that e-05 renders that value close to 0.

cbind(rownames(rf_model$fit$importance), rf_model$fit$importance) %>% as_tibble() %>% select(predictor= V1, MeanDecreaseAccuracy)   %>% 
  mutate(MeanDecreaseAccuracy=as.numeric(MeanDecreaseAccuracy)) %>% arrange(desc((MeanDecreaseAccuracy))) %>% head()

# As the model is trained using the randomForest package, we can use randomForest::importance to extract the results too. 
# Use argument type=1 to extract permutation-based importance and set scale=F to prevent normalization. 
# There are benefits for not scaling the permutation-based score.

cbind(rownames(randomForest::importance(rf_model$fit, type=1, scale=F)) %>% as_tibble(), randomForest::importance(rf_model$fit, type=1, scale=F) %>% as_tibble())%>% arrange(desc(MeanDecreaseAccuracy)) %>% head()

# importance of each predictor
important <- randomForest::importance(rf_model$fit, type=1, scale=F)
Important_Features <- data.frame(Feature = row.names(important), Importance = important [,1]) #to plot the Mean Decrease Accuracy
write.csv(Important_Features, "Important_Features_model4.csv")
Important_Features <- read.csv("Important_Features_model4.csv")

#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               Importance) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy")+
  ggtitle("Important Features in Random Forest\n") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 6))
plot_

# 2.2- Impurity-based importance

# Likewise, we can use both approaches to calculate impurity-based importance. Using the tidyverse way:
cbind(rownames(rf_model$fit$importance), rf_model$fit$importance) %>% as_tibble() %>% select(predictor= V1, MeanDecreaseGini)   %>% 
  mutate(MeanDecreaseGini=as.numeric(MeanDecreaseGini)) %>% arrange(desc((MeanDecreaseGini))) %>% head(10)

# For the randomForest::importance function, we will change the argument type to 2. The scale argument does not work for impurity-based scores.
cbind(rownames(randomForest::importance(rf_model$fit, type=2)) %>% as_tibble(), randomForest::importance(rf_model$fit, type=2) %>% as_tibble())%>% arrange(desc(MeanDecreaseGini)) %>% head(10)

# importance of each predictor
gini <- randomForest::importance(rf_model$fit, type=2, scale=F)
Gini_Features <- data.frame(Feature = row.names(gini), Importance = gini [,1]) #to plot the Mean Decrease Accuracy
write.csv(Gini_Features, "Gini_Features_model4.csv")
Gini_Features <- read.csv("Gini_Features_model4.csv")

#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Gini_Features, 
                aes(x= reorder(Feature,
                               Importance) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Gini") +
  ggtitle("Important Features in Random Forest\n") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 6))
plot_

# In this post we did a post-hoc analysis on random forest to understand the model by using permutation and impurity variable importance ranking. 
# In the next post, we will continue post-hoc analysis on random forest model using other techniques.


# Reference tutorial: https://www.r-bloggers.com/2019/08/explaining-predictions-random-forest-post-hoc-analysis-randomforestexplainer-package/
# Explaining Predictions: Random Forest Post-hoc Analysis (randomForestExplainer package)

#1- Measuring variable importance

# The randomForestExplainer generates more feature importance score than the randomForest package. 
# These scores generated by the randomForestExplainer::measure_importance function can be categorized into 3 types of feature importance score.

# 1.1- Measures of importance based on the structure of the forest which is not examine in the randomForest package.
# i) Mean_minimal_depth: Minimal depth for a variable in a tree equals to the depth of the node which splits on that variable and is the closest to the root of the tree. 
# If it is low then a lot of observations are divided into groups on the basis of this variable.
# ii) no_of_nodes: Usually, number of trees and number of nodes are the same if trees are shallow. The number of nodes contain similar information as the number of trees and p-values. Thus, no_of_tress and p_value are omitted in plot_importance_ggpairs graph (which will be described later).
# iii) no_of_trees
# iv) p_value: Number of nodes in which predictor X was used for splitting exceeds the theoretical number of successes if they were random, following the binomial distribution given.
# v) times_a_root: Measures of importance based on the decrease in predictive accuracy post variable perturbation

# 1.2- 
# i) accuracy_decrease: Use for classification problems. Same value as randomForest::importance(rf_model$fit, type=2)
# ii) mse_increase: Use for regression problems.

# 1.3- Measure of importance based on loss function
# i) gini_decrease Use for classifcation cases. Same value as randomForest::importance(rf_model$fit, type=1).
#ii) node_purity_increase: Use for regression cases. Measured by the decrease in sum of squares

library(randomForestExplainer)
impt_frame<-measure_importance(rf_model$fit)
impt_frame %>% head()

# 1.4- Plotting
# plot_multi_way_importance plots variable importance scores from importance_frame. 
# By default, the function plots scores measuring importance based on the structure of the forest, mean_min_depth against times_a_root. 
plot_multi_way_importance(impt_frame, no_of_labels = 6)

# plot_multi_way_importance can plot up to 3 feature importance scores using any of the three kinds of importance measurements.

plot_multi_way_importance(impt_frame, x_measure = "accuracy_decrease", y_measure = "gini_decrease", size_measure = "p_value", no_of_labels = 9)

# 1.5- Selecting optimal importance scores to plot

# There are many multi-way importance plots which one can create thus we need to identify which sub-set of importance scores will provide us with more helpful plots. 
# Examining the relationship between the importance measurements with plot_importance_ggpairs can assist us in identifying a sub-set of scores to use. 
# We can pick three scores that least agree with each other, points in plots which are most dispersed.

plot_importance_ggpairs(impt_frame)
plot_importance_rankings(impt_frame)

# 1.6- Variable depth
# Previously, we looked at various types of importance measures. 
# Now we are specifically examining the mean minimal depth in detail. 
# The distribution of the mean minimal depth allows us to appreciate the variable’s role in the random forest’s structure and prediction. 
# plot_min_depth_distribution plots the top ten variables according to mean minimal depth calculated using top trees. 
# The mean minimal depth can be calculated in 3 different ways in plot_min_depth_distribution using the mean_sample argument. 
# The calculation differs in the way they treat missing values that appear when a feature is not used for tree splitting. 
# As a result, the ranking of variables may change for each calculation.

# The mean minimal depth is indicated by a vertical bar with the mean value beside it. 
# The smaller the mean minimal depth, the more important the variable is and the higher up the y-axis the variable will be. 
# The rainbow gradient reveals the min and max minimal depth for each variable. 
# The bigger the proportion of minimal depth zero (red blocks), the more frequent the variable is the root of a tree. 
# The smaller the proportion of NA minimal depth (gray blocks), the more frequent the variable is used for splitting trees. 
# The range of the x-axis is from zero to the maximum number of trees for the feature.
md_frame <- min_depth_distribution(rf_model$fit)
plot_min_depth_distribution(md_frame, mean_sample = "top_trees") # default mean_sample arg 

# 1.7- Variable Interaction
# We will use the important_variables function to select the top 6 variables based on the following variable importance measurements, times_a_root and no_of_nodes.
vars<- important_variables(impt_frame, k = 6, measures = c("times_a_root", "no_of_nodes"))

# After identifying the top 6 variables, we can examine the interactions between the variables with the min_depth_interactions function.
# The interaction is reflected as the mean_min_depth which is the mean conditional minimal depth, where a variable is taken as a root node/root_variable and the mean minimal depth is calculated for the other variable. 
# The uncond_mean_min_depth represents the unconditional mean minimal depth of the variable which is the mean minimal depth of the variable without having a stipulated root variable. 
# This value is the same as the mean value seen on the vertical bar in plot_min_depth_distribution.
interactions_frame <- min_depth_interactions(rf_model$fit, vars)
plot_min_depth_interactions(interactions_frame)

# 1.8- Interactive variables and forest prediction

# We can further evaluate the variable interactions by plotting the probability of a prediction against the variables making up the interaction. 
# For instance we plot the probability of having heart disease against resting blood pressure rest_bp and ST depression duration during exercise test ex_STdepression_dur. 
# The interaction of these two variables are the most frequent interaction as seen in plot_min_depth_interactions. 
# We plot the forest prediction against interactive variables with plot_predict_interaction.

plot_predict_interaction(rf_model$fit, bake(covid_prep, new_data = train_data), "IL12p70", "PD_L1")

set.seed(1)
forest <- randomForest::randomForest(Outcome ~ ., data = bake(covid_prep, new_data = train_data), localImp = TRUE, importance=T, ntree = 1000, type= "classification")

predict_plot<-plot_predict_interaction(forest, bake(covid_prep, new_data = train_data), "IL12p70", "PD_L1", main = "Distribution of Predicted Probability of Fatal COVID-19 Outcome") + 
  theme(legend.position="bottom") + 
  geom_hline(yintercept = 400, linetype="longdash") + 
  geom_vline(xintercept = 40, linetype="longdash")

test_plot<-ggplot(test_data, aes(IL12p70, PD_L1, colour=Outcome)) +geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", y="PD-L1") + 
  theme_light(base_size = 16) +  
  scale_color_manual(values=c("dark grey", "dark blue")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 400, linetype="dotted") +
  geom_vline(xintercept = 40, linetype="dotted")
test_plot

gridExtra::grid.arrange(predict_plot, test_plot, ncol=2)

predict_plot<-plot_predict_interaction(forest, bake(covid_prep, new_data = train_data), "IL12p70", "Creatinine", main = "Distribution of Predicted Probability of Fatal COVID-19 Outcome") + 
  theme(legend.position="bottom") + 
  geom_hline(yintercept = 1.5, linetype="longdash") + 
  geom_vline(xintercept = 40, linetype="longdash")

test_plot<-ggplot(test_data, aes(IL12p70, Creatinine, colour=Outcome)) +geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", y="Creatinine") + 
  theme_light(base_size = 16) + 
  scale_color_manual(values=c("dark grey", "dark blue")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 1.5, linetype="dotted") +
  geom_vline(xintercept = 40, linetype="dotted")
test_plot

gridExtra::grid.arrange(predict_plot, test_plot, ncol=2)

# 2- Run random forest rmodel with the variables important for prediction as determined above
# 2.1- Recovered vs fatal
# train model with rand_forest
set.seed(1)
rf_model_final<-rand_forest(trees = 1000,  mode = "classification") %>% 
  set_engine("randomForest",  importance=T, localImp = T, ) %>% # importance = T to have permutation score calculated
  fit(Outcome ~ IL12p70+PD_L1+IFNb+IFNlambda2+IL33+IL27+CCL1+Platelets+Creatinine, data = juice(covid_prep)) # localImp=T for randomForestExplainer(next post)
saveRDS(rf_model_final,"rftrain_recovered_fatal.rds")
rf_model_final <- readRDS("rftrain_recovered_fatal.rds")

plot(rf_model_final)
reprtree:::plot.getTree(rf_model_final, k=1)
# summarize the fit
summary(rf_model_final)
# make predictions in the test data set
predictions <- predict(rf_model_final, test_data)
# summarize accuracy and predicitons in the test data set
table(predictions, test_data$Outcome)
mean(predictions == test_data$Outcome)
confusionMatrix(predictions, test_data$Outcome)

# train model with random_forest
set.seed(7)
forest_final <- randomForest::randomForest(Outcome ~ IL12p70+PD_L1+IFNb+IFNlambda2+IL33+IL27+CCL1+Platelets+Creatinine, data = bake(covid_prep, new_data = train_data), localImp = TRUE, importance=T, ntree = 1000, type= "classification")

saveRDS(forest_final,"forest_final_recovered_fatal.rds")
forest_final <- readRDS("forest_final_recovered_fatal.rds")

plot(forest_final)
reprtree:::plot.getTree(forest_final, k=1)
# summarize the fit
summary(forest_final)
# make predictions in the test data set
predictions <- predict(forest_final, test_data)
# summarize accuracy and predicitons in the test data set
table(predictions, test_data$Outcome)
mean(predictions == test_data$Outcome)
confusionMatrix(predictions, test_data$Outcome)

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



