library(rpart)
library(readr)
library(caTools)
library(dplyr)
library(party)
library(partykit)
library(rpart.plot)
library(C50)
library(randomForest)
library(reprtree)
library(e1071)
library(caret)
library(htree)
library(LongituRF)
library(VSURF)
library(tidyverse)
library(tidymodels)
library(pROC)

# reading file (importar o arquivo)
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_htree.csv")

# data frame para cluster
mydata <- as.data.frame(multiplex)
mydata$Outcome <- as.factor(mydata$Outcome)


mydata2 <- mydata[1:47,] #RF analysis within the fatal group for prediction of disease progression
mydata2$Outcome2 <- as.factor(mydata2$Outcome2)

# 1- RF longitudinal data with the htree package

# fit historical random forest 
# (NOTE: response is 0/1 so a regression tree is the same as a classification tree with Gini-index splitting)
# -- set concurrent and historical predictors
historical_predictors=match(c("SCF", "OPN",	"Ang1",	"Ang2",	"ANG2_ANG1",	"ICAM1",	"VCAM1",	"E_selectin",
                              "Syndecan1",	"ADAMTS13",	"vWFA2",	"TF",	"Fibronectin",	"D_Dimer",	"CD40L",	"CXCL4",	"CXCL7",	"TPO",	
                              "IL11",	"IL1a",	"IL1b",	"IL18",	"IL33",	"IL1RA",	"TNFa", "IL6",	"IL8",	
                              "IL12p70",	"IL17",	"IL4",	"IL5",	"IL10", "IL13",	"IL7",	
                              "IL15",	"IL21",	"CCL11",	"CCL24",	"G_CSF",	"GM_CSF",
                              "M_CSF",	"IL27",	"MPO",	"L_selectin",	"CCL1",	"MCP1",	"MCP2",	"MIP1a",	"MIP1b",
                              "CXCL1",	"CXCL2",	"MIP3a",	"MIP3b",	"RANTES",	"CXCL9",	"CXCL10",	"CXCL11",	"CXC3CL1",	"IFNg",	"IFNa",
                              "IFNb",	"IFNlambda2",	"IFNlambda3",	"PD_L1",	"FasL",	"GranzymeB",	"EGF",	"FGFb",	"HGF",	"PDGF_AA",	"PDGF_BB",
                              "C9",	"Hemoglobin",	"Hematocrit",	"Lymphocytes", "Neutrophils", 	"NCLR",	"Platelets",	"ALT",	"Glucose",	"Creatinine",	
                              "Urea",	"LDH",	"CRP"), names(mydata))

control=list(vh=historical_predictors, ntrees=1000, nodesize=5, classify=TRUE, mtry=5)

htree = hrf(x=mydata, id=mydata$RecordID, time=mydata$Day, yindx="Outcome", control=control, se=TRUE) # model for the outcome recovered vs fatal
htree2 = hrf(x=mydata2, id=mydata2$RecordID, time=mydata2$Day, yindx="Outcome2", control=control) # model for early vs late death

# out-of-bag error
plot(1:length(htree$error),htree$error,type="l",main="OOB error",xlab="forest size",ylab="mse")
points(1:length(htree$error),htree$error,type="l",col="blue")

# -- variable importance table 
set.seed(123)
vi=varimp_hrf(htree) 
vi
select_row <- rownames(vi)
select_row <- select_row[!select_row %in% c("1","86", "87", "88","89", "90")]
vi <- vi[select_row,]

Important_Features <- data.frame(Feature = vi$Predictor, Importance = vi$`Z-value`)
#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = (Importance))) +
  geom_bar(stat = "identity", 
           fill = "#800080") +
  coord_flip() +
  theme_light(base_size = 8) +
  xlab("") + 
  ylab("Importance") +
  ggtitle("Important Features in Random Forest\n") +
  theme(plot.title = element_text(size=10))
plot_

vi
write.csv(vi, "COVID_Clinical_Luminex_deceased_recovered_longitudinal_htree_VI_list.csv")

p=predict_hrf(htree, time=mydata$Day, id=mydata$RecordID)
p <- cbind(p, mydata$RecordID)
write.csv(p, "COVID_Clinical_Luminex_deceased_recovered_longitudinal_htree_predictions.csv")

sample_data = sample.split(mydata, SplitRatio = 0.7)
train_data <- subset(mydata, sample_data == TRUE)
test_data <- subset(mydata, sample_data == FALSE)
mydata2 <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_htree.csv")
p2=predict_hrf(htree, test_data, time=test_data$Day, id=test_data$RecordID)
p2 <- cbind(p2, test_data$RecordID)
write.csv(p3, "COVID_Clinical_Luminex_deceased_recovered_longitudinal_htree_predictions_testdata.csv")

# Analysis for the early vs late death prediction
# out-of-bag error
plot(1:length(htree2$error),htree2$error,type="l",main="OOB error",xlab="forest size",ylab="mse")
points(1:length(htree2$error),htree2$error,type="l",col="blue")

# -- variable importance table 
set.seed(123)
vi=varimp_hrf(htree2) 
vi
select_row <- rownames(vi)
select_row <- select_row[!select_row %in% c("1","86", "87", "88","89", "90")]
vi <- vi[select_row,]

Important_Features <- data.frame(Feature = vi$Predictor, Importance = vi$`Z-value`)
#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = (Importance))) +
  geom_bar(stat = "identity", 
           fill = "#800080") +
  coord_flip() +
  theme_light(base_size = 8) +
  xlab("") + 
  ylab("Importance") +
  ggtitle("Important Features in Random Forest\n") +
  theme(plot.title = element_text(size=10))
plot_

p=predict_hrf(htree, time=mydata$Day, id=mydata$RecordID)
p <- cbind(p, mydata$RecordID)
write.csv(p, "COVID_Clinical_Luminex_deceased_recovered_longitudinal_htree_predictions.csv")

sample_data = sample.split(mydata, SplitRatio = 0.7)
train_data <- subset(mydata, sample_data == TRUE)
test_data <- subset(mydata, sample_data == FALSE)
mydata2 <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_htree.csv")
p2=predict_hrf(htree, test_data, time=test_data$Day, id=test_data$RecordID)
p2 <- cbind(p2, test_data$RecordID)
write.csv(p3, "COVID_Clinical_Luminex_deceased_recovered_longitudinal_htree_predictions_testdata.csv")

#####################################################################################################################################################################################################################################################################################################################################################################

library(LongituRF)

# 2- RF longitudinal data with the LongituRF package - Recovered vs Fatal
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_vsurf.csv")

mydata <- multiplex

# Split the train and test datasets
set.seed(123)
sample_data <- initial_split(mydata, prop=0.7, strata = "Outcome2")
train_data <- training(sample_data)
test_data <- testing(sample_data)

# Trying as in https://www.christopherloan.com/blog/running-longitudinal-random-forests-with-longiturf/

# 2.1- Organizing the data
list_dat <- 
  list(
    # predictors, level1 and level2
    X = train_data[,2:85], 
    # outcome
    Y = as.numeric(pull(train_data, Outcome2)), 
    # id variables for each unique individual (?)
    id = as.numeric(pull(train_data, RecordID)), 
    # random effects (using only a random intercept)
    Z = as.matrix(rep(1, nrow(train_data))), 
    # years where wave 1 = 0, wave 2 = 1, wave 3 = 2
    time = as.numeric(pull(train_data, Day)) 
  )

# See example of the data structure below
data <- DataLongGenerator(n=17, p=6,G=6)

list_dat <- 
  list(
    # predictors, level1 and level2
    X = mydata[,2:85], 
    # outcome
    Y = as.numeric(pull(mydata, Outcome2)), 
    # id variables for each unique individual (?)
    id = as.numeric(pull(mydata, RecordID)), 
    # random effects (using only a random intercept)
    Z = as.matrix(rep(1, nrow(mydata))), 
    # years where wave 1 = 0, wave 2 = 1, wave 3 = 2
    time = as.numeric(pull(mydata, Day)) 
  )


# 2.2- Training data with the different algorithms

# MERT
start_time <- proc.time()
mert1 <- MERT(X = data.frame(list_dat$X), Y = list_dat$Y,
    id = list_dat$id, Z = list_dat$Z, time = list_dat$time, sto = 'OrnUhl') # Ornstein–Uhlenbeck
    #iter = 100, 
    #delta = 0.001

# MERF
start_time <- proc.time()
merf1 <- MERF(X = data.frame(list_dat$X), Y = list_dat$Y, id = list_dat$id, Z = list_dat$Z,
    time = list_dat$time, sto = 'none')
    #iter = 100,
    #delta = 0.001
    #mtry = ceiling(ncol(data.frame(list_dat$X))/3),
    #ntree = 500
    
# REEMTree
start_time <- proc.time()
reemtree1 <- REEMtree(X = data.frame(list_dat$X), Y = list_dat$Y, id = list_dat$id,
    Z = list_dat$Z, time = list_dat$time, sto = 'OrnUhl')
    
# REEMforest  
# https://stackoverflow.com/questions/66188123/implementing-longitudinal-random-forest-with-longiturf-package-in-r it never works!
start_time <- proc.time()
reemforest1 <- REEMforest(X = data.frame(list_dat$X), Y = list_dat$Y, id = list_dat$id, mtry = 3,
                      Z = list_dat$Z, time = list_dat$time, sto = "OrnUhl")    
    
# 2.3- Variable Importance for feature selection - based on trained data

# MERT
tibble(
  variables = names(mert1$forest$variable.importance), 
  importance = mert1$forest$variable.importance
) %>% 
  ggplot() + geom_col(aes(y = fct_reorder(variables, importance), x = importance),
    fill = 'coral3', alpha = 0.5, color = 'black', show.legend = FALSE) + labs(title = 'MERT', y = element_blank()) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

Important_Features_MERT <- data.frame(mert1$forest$variable.importance)
write.csv(Important_Features_MERT, "Important_Features_MERT_train_data.csv")

Important_Features_MERT <- read.csv("Important_Features_MERT_train_data.csv")
# use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features_MERT, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) + geom_bar(stat = "identity") +
  coord_flip() + theme_bw(base_size = 10) +
  xlab("") +  ylab("Ikportance") + ggtitle("MERT") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 10))
plot_

# REEMTree

tibble(
  variables = names(reemtree1$forest$variable.importance), 
  importance = reemtree1$forest$variable.importance
) %>% 
  ggplot() + 
  geom_col(aes(y = fct_reorder(variables, importance), x = importance), fill = 'orange', alpha = 0.5, color = 'black', show.legend = FALSE) +
  labs(title = 'REEMtree', y = element_blank()) +
  theme_minimal() + theme(axis.text.x = element_blank())

# MERF
tibble(
  variables = names(merf1$forest$importance[,1]), 
  importance = merf1$forest$importance[,1]
) %>% 
  ggplot() + 
  geom_col(aes(y = fct_reorder(variables, importance), x = importance), fill = 'darkgreen', alpha = 0.5, color = 'black', show.legend = FALSE) +
  labs(title = 'MERF', y = element_blank()) +
  theme_minimal() + theme(axis.text.x = element_blank())

Important_Features_MERF <- data.frame(merf1$forest$importance[,1])
write.csv(Important_Features_MERF, "Important_Features_MERF_train_data.csv")

Important_Features_MERF <- read.csv("Important_Features_MERF_train_data.csv")
# use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features_MERF, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 10) +
  xlab("") + 
  ylab("Ikportance") +
  ggtitle("MERF") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 8))
plot_

# 2.4- Making Predictions on Unseen Data

list_test <- 
  list(
    # predictors, level1 and level2
    X = test_data[,2:85], 
    # outcome
    Y = as.numeric(pull(test_data, Outcome2)), 
    # id variables for each unique individual (?)
    id = as.numeric(pull(test_data, RecordID)), 
    # random effects (using only a random intercept)
    Z = as.matrix(rep(1, nrow(test_data))),
    # years where wave 1 = 0, wave 2 = 1, wave 3 = 2
    time = as.numeric(pull(test_data, Day)) 
  )

# Obtain predictions for each of the three LongituRF based models (the (S)MERT, (S)REEMtree, and the (S)MERF). 
# The native prediction function from LongituRF removed the predictions for cases that weren’t observed in the data.
fixed_prediction_function <- function(object, X, id, Z, time, new_df, id_var_name, ...) {
  `%notin%` <- Negate(`%in%`)
    if("tidyverse" %notin% (.packages())){suppressMessages(library(tidyverse))}
    
    preds_existing <- 
      predict(
        object = object,
        X = X,
        id = id,
        Z = Z,
        time = time
      )
    
    temp <- 
      new_df %>% 
      filter({{id_var_name}} %notin% object$id) %>% 
      mutate(predictions = predict(object = object$forest, newdata = .))
    
    final_df <- 
      new_df %>% 
      filter({{id_var_name}} %in% object$id) %>% 
      mutate(predictions = preds_existing) %>% 
      bind_rows(temp)
    return(final_df)
}

predictions_mert_df <- 
  fixed_prediction_function(
    object = reemtree1,
    X = data.frame(list_test$X),
    id = list_test$id,
    Z = list_test$Z,
    time = list_test$time)

############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

library(LongituRF)

# 3- RF longitudinal data with the LongituRF package - ED vs LD
multiplex <- read.csv("COVID_Clinical_Luminex_ED_LD_simon_missing_data_imputed.csv")
mydata <- multiplex

# Split the train and test datasets
set.seed(123)
sample_data <- initial_split(mydata, prop=0.7, strata = "Outcome")
train_data <- training(sample_data)
test_data <- testing(sample_data)

# Trying as in https://www.christopherloan.com/blog/running-longitudinal-random-forests-with-longiturf/

# 3.1- Organizing the data
list_dat <- 
  list(
    # predictors, level1 and level2
    X = mydata[,2:85], 
    # outcome
    Y = as.numeric(pull(mydata, Outcome)), 
    # id variables for each unique individual (?)
    id = as.numeric(pull(mydata, RecordID)), 
    # random effects (using only a random intercept)
    Z = as.matrix(rep(1, nrow(mydata))), 
    # years where wave 1 = 0, wave 2 = 1, wave 3 = 2
    time = as.numeric(pull(mydata, Day)) 
  )

# ignore below
##list_dat <- 
 # list(
    # predictors, level1 and level2
   # X = train_data[,2:85], 
    # outcome
  #  Y = as.numeric(pull(train_data, Outcome)), 
    # id variables for each unique individual (?)
  #  id = as.numeric(pull(train_data, RecordID)), 
    # random effects (using only a random intercept)
  #  Z = as.matrix(rep(1, nrow(train_data))), 
    # years where wave 1 = 0, wave 2 = 1, wave 3 = 2
  #  time = as.numeric(pull(train_data, Day)) 
 # )


# 3.2- Training data with the different algorithms
# MERT
start_time <- proc.time()
mert1 <- MERT(X = data.frame(list_dat$X), Y = list_dat$Y,
              id = list_dat$id, Z = list_dat$Z, time = list_dat$time, sto = 'OrnUhl') # Ornstein–Uhlenbeck
#iter = 100, 
#delta = 0.001

# MERF
start_time <- proc.time()
merf1 <- MERF(X = data.frame(list_dat$X), Y = list_dat$Y, id = list_dat$id, Z = list_dat$Z,
              time = list_dat$time, sto = 'none')
#iter = 100,
#delta = 0.001
#mtry = ceiling(ncol(data.frame(list_dat$X))/3),
#ntree = 500

# REEMTree
start_time <- proc.time()
reemtree1 <- REEMtree(X = data.frame(list_dat$X), Y = list_dat$Y, id = list_dat$id,
                      Z = list_dat$Z, time = list_dat$time, sto = 'OrnUhl')

# 3.3- Variable Importance for feature selection - based on trained data

# MERT
tibble(
  variables = names(mert1$forest$variable.importance), 
  importance = mert1$forest$variable.importance
) %>% 
  ggplot() + geom_col(aes(y = fct_reorder(variables, importance), x = importance),
                      fill = 'coral3', alpha = 0.5, color = 'black', show.legend = FALSE) + labs(title = 'MERT', y = element_blank()) +
  theme_minimal() +
  theme(axis.text.x = element_blank())

Important_Features_MERT <- data.frame(mert1$forest$variable.importance)
write.csv(Important_Features_MERT, "Important_Features_MERT_train_data.csv")

Important_Features_MERT <- read.csv("Important_Features_MERT_train_data.csv")
# use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features_MERT, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) + geom_bar(stat = "identity") +
  coord_flip() + theme_bw(base_size = 10) +
  xlab("") +  ylab("Ikportance") + ggtitle("MERT") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 10))
plot_

# REEMTree
tibble(
  variables = names(reemtree1$forest$variable.importance), 
  importance = reemtree1$forest$variable.importance
) %>% 
  ggplot() + 
  geom_col(aes(y = fct_reorder(variables, importance), x = importance), fill = 'orange', alpha = 0.5, color = 'black', show.legend = FALSE) +
  labs(title = 'REEMtree', y = element_blank()) +
  theme_minimal() + theme(axis.text.x = element_blank())

# MERF
tibble(
  variables = names(merf1$forest$importance[,1]), 
  importance = merf1$forest$importance[,1]
) %>% 
  ggplot() + 
  geom_col(aes(y = fct_reorder(variables, importance), x = importance), fill = 'darkgreen', alpha = 0.5, color = 'black', show.legend = FALSE) +
  labs(title = 'MERF', y = element_blank()) +
  theme_minimal() + theme(axis.text.x = element_blank())

Important_Features_MERF <- data.frame(merf1$forest$importance[,1])
write.csv(Important_Features_MERF, "Important_Features_MERF_train_data.csv")

Important_Features_MERF <- read.csv("Important_Features_MERF_train_data.csv")
# use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features_MERF, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 10) +
  xlab("") + 
  ylab("Ikportance") +
  ggtitle("MERF") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 8))
plot_


#####################################################################################################################################################################################################################################################################################################################################################################
# 4- Variable Importance and Selection with VSURF package
# 4.1- Recovered vs Fatal

multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_vsurf.csv")
mydata <- multiplex[,-1]
mydata2 <- mydata[,1:85]
mydata2 <- as.data.frame(mydata2)

set.seed(1)
vsurf <- VSURF(Outcome~., data=mydata2, mtry=5, ntree=1000)
summary(vsurf)

plot(vsurf,  step = "all", var.names = TRUE,
     imp.mean = TRUE, imp.sd = TRUE,
     nvar.imp.mean = length(vsurf$imp.mean.dec),
     nvar.imp.sd = length(vsurf$imp.sd.dec),
     nvar.interp = length(vsurf$varselect.thres),
     nvar.pred = length(vsurf$varselect.pred))

# The computation of the 50 forests, the ranking and elimination steps are obtained with the VSURF_thres function:
set.seed(123)
vsurf.thres <- VSURF_thres(Outcome~., data=mydata2, mtry = 5, ntree=1000)
vsurf.thres$varselect.thres

plot(vsurf, step = "thres", imp.mean = FALSE, var.names=TRUE)
plot(vsurf.thres, var.names = TRUE, imp.mean = TRUE,
     imp.sd = TRUE, nvar.imp.mean = length(vsurf$imp.mean.dec),
     nvar.imp.sd = length(vsurf$imp.sd.dec))

Important_Features <- data.frame(Feature = vsurf.thres$varselect.thres, Importance = vsurf.thres$imp.varselect.thres)
write.csv(Important_Features, "Important_Features_VSURF2.csv")

ortho <- read.csv("Important_Features_VSURF2.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, file="Important_Features_VSURF2.csv")

Important_Features <- read.csv("Important_Features_VSURF_final.csv")
# use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 10) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for thresholding") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 10))
plot_

#Variable selection procedure for interpretation.
set.seed(123)
vsurf.interp <- VSURF_interp(Outcome~., data=mydata2, vars = vsurf.thres$varselect.thres)
vsurf.interp$varselect.interp

#method for class 'VSURF_interp'
plot(vsurf.interp, var.names = TRUE,
     nvar.interp = length(vsurf.interp$varselect.interp))

#use ggplot2 to plot our set of important features for interpretation
vsurf.interp2 <- vsurf.interp$err.interp[1:8]
Important_Features <- data.frame(Feature = vsurf.interp$varselect.interp, Importance = vsurf.interp2)
write.csv(Important_Features, "Interpretation_Features_VSURF2.csv")

ortho <- read.csv("Interpretation_Features_VSURF2.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Interpretation_Features_VSURF2.csv")

Important_Features <- read.csv("Interpretation_Features_VSURF2.csv")
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for interpretation") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 16))
plot_

#Variable selection procedure for prediction.
set.seed(123)
vsurf.pred <- VSURF_pred(Outcome~., data=mydata2, err.interp = vsurf.interp$err.interp, varselect.interp = vsurf.interp$varselect.interp)
plot(vsurf.pred, var.names = TRUE,
     nvar.pred = length(vsurf.pred$varselect.pred))

#use ggplot2 to plot our set of important features for prediction:
Important_Features <- data.frame(Feature = vsurf.pred$varselect.pred, Importance = vsurf.pred$err.pred)
write.csv(Important_Features, "Prediction_Features_VSURF2.csv")

ortho <- read.csv("Prediction_Features_VSURF2.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Prediction_Features_VSURF2.csv")

Important_Features <- read.csv("Interpretation_Features_VSURF_final.csv")
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for prediction") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 16))
plot_

#####################################################################################################################################################################################################################################################################################################################################################################

# 4.2- Variable Importance and Selection with VSURF package - Recovered vs Early vs Late

multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_all_samples_vsurf.csv")
mydata <- multiplex[,-1]
mydata2 <- mydata[,1:85]
mydata2 <- as.data.frame(mydata2)

set.seed(123)
vsurf <- VSURF(Outcome~., data=mydata2, mtry=5, ntree=1000)
summary(vsurf)
plot(vsurf,  step = "all", var.names = TRUE,
     imp.mean = TRUE, imp.sd = TRUE,
     nvar.imp.mean = length(vsurf$imp.mean.dec),
     nvar.imp.sd = length(vsurf$imp.sd.dec),
     nvar.interp = length(vsurf$varselect.thres),
     nvar.pred = length(vsurf$varselect.pred))

#The computation of the 50 forests, the ranking and elimination steps are obtained with the VSURF_thres function:
vsurf.thres <- VSURF_thres(Outcome~., data=mydata2, mtry = 5, ntree=1000)
vsurf.thres$varselect.thres

plot(vsurf, step = "thres", imp.mean = FALSE, var.names=TRUE)
plot(vsurf.thres, var.names = TRUE, imp.mean = TRUE,
     imp.sd = TRUE, nvar.imp.mean = length(vsurf$imp.mean.dec),
     nvar.imp.sd = length(vsurf$imp.sd.dec))

Important_Features <- data.frame(Feature = vsurf.thres$varselect.thres, Importance = vsurf.thres$imp.varselect.thres)
write.csv(Important_Features, "Important_Features_Recov_ED_LD.csv")

ortho <- read.csv("Important_Features_Recov_ED_LD.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Important_Features_Recov_ED_LD.csv")

Important_Features <- read.csv("Important_Features_Recov_ED_LD.csv")
#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light(base_size = 10) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for thresholding") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 10))
plot_

#Variable selection procedure for interpretation.
vsurf.interp <- VSURF_interp(Outcome~., data=mydata2, vars = vsurf.thres$varselect.thres)
vsurf.interp$varselect.interp

#method for class 'VSURF_interp'
plot(vsurf.interp, var.names = TRUE,
     nvar.interp = length(vsurf.interp$varselect.interp))

#use ggplot2 to plot our set of important features for interpretation
vsurf.interp2 <- vsurf.interp$err.interp[1:8]
Important_Features <- data.frame(Feature = vsurf.interp$varselect.interp, Importance = vsurf.interp2)
write.csv(Important_Features, "Interpretation_Features_Recov_ED_LD.csv")

ortho <- read.csv("Interpretation_Features_Recov_ED_LD.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Interpretation_Features_Recov_ED_LD.csv")

Important_Features <- read.csv("Interpretation_Features_Recov_ED_LD.csv")
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for interpretation") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 16))
plot_

#Variable selection procedure for prediction.
vsurf.pred <- VSURF_pred(Outcome~., data=mydata2, err.interp = vsurf.interp$err.interp, varselect.interp = vsurf.interp$varselect.interp)
plot(vsurf.pred, var.names = TRUE,
     nvar.pred = length(vsurf.pred$varselect.pred))

#use ggplot2 to plot our set of important features for prediction:
Important_Features <- data.frame(Feature = vsurf.pred$varselect.pred, Importance = vsurf.pred$err.pred)
write.csv(Important_Features, "Prediction_Features_Recov_ED_LD.csv")

ortho <- read.csv("Prediction_Features_Recov_ED_LD.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Prediction_Features_Recov_ED_LD.csv")

Important_Features <- read.csv("Prediction_Features_Recov_ED_LD.csv")
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for prediction") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 16))
plot_

############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 4.3- Variable Importance and Selection with VSURF package - Early vs Late Death only 
# tested with or without imputed missing data

multiplex <- read.csv("COVID_Clinical_Luminex_ED_LD_simon_missing_data.csv")
mydata2 <- multiplex[,-1]
#mydata2 <- mydata[,1:85]
#mydata2 <- as.data.frame(mydata2)
#mydata2 <- mydata2[1:47,]

#missing data in RF: https://stackoverflow.com/questions/8370455/how-to-use-random-forests-in-r-with-missing-values
set.seed(222)
mydata2.imputed <- rfImpute(Outcome ~ ., mydata2, iter=10, ntree=1000)
write.csv(mydata2.imputed, "COVID_Clinical_Luminex_ED_LD_simon_missing_data_imputed.csv")
mydata2 <- read.csv("COVID_Clinical_Luminex_ED_LD_simon_missing_data_imputed.csv")
mydata2 <- mydata2[,-1]

set.seed(1234)
vsurf <- VSURF(Outcome~., data=mydata2, mtry = 5, ntree=1000)
summary(vsurf)
plot(vsurf,  step = "all", var.names = TRUE,
     imp.mean = TRUE, imp.sd = TRUE,
     nvar.imp.mean = length(vsurf$imp.mean.dec),
     nvar.imp.sd = length(vsurf$imp.sd.dec),
     nvar.interp = length(vsurf$varselect.thres),
     nvar.pred = length(vsurf$varselect.pred))

#The computation of the 50 forests, the ranking and elimination steps are obtained with the VSURF_thres function:
set.seed(123)
vsurf.thres <- VSURF_thres(Outcome~., data=mydata2, mtry = 5, ntree=1000)
vsurf.thres$varselect.thres

plot(vsurf, step = "thres", imp.mean = FALSE, var.names=TRUE)
plot(vsurf.thres, var.names = TRUE, imp.mean = TRUE,
     imp.sd = TRUE, nvar.imp.mean = length(vsurf$imp.mean.dec),
     nvar.imp.sd = length(vsurf$imp.sd.dec))

Important_Features <- data.frame(Feature = vsurf.thres$varselect.thres, Importance = vsurf.thres$imp.varselect.thres)
write.csv(Important_Features, "Important_Features_ED_LD.csv")

ortho <- read.csv("Important_Features_ED_LD.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Important_Features_ED_LD.csv")

Important_Features <- read.csv("Important_Features_ED_LD.csv")
#use ggplot2 to plot our set of important features:
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 10) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for thresholding") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 10))
plot_

#Variable selection procedure for interpretation.
vsurf.interp <- VSURF_interp(Outcome~., data=mydata2, vars = vsurf.thres$varselect.thres)
vsurf.interp$varselect.interp

#method for class 'VSURF_interp'
plot(vsurf.interp, var.names = TRUE,
     nvar.interp = length(vsurf.interp$varselect.interp))

#use ggplot2 to plot our set of important features for interpretation
vsurf.interp2 <- vsurf.interp$err.interp[1:5]
Important_Features <- data.frame(Feature = vsurf.interp$varselect.interp, Importance = vsurf.interp2)
write.csv(Important_Features, "Interpretation_Features_ED_LD.csv")

ortho <- read.csv("Interpretation_Features_ED_LD.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Interpretation_Features_ED_LD.csv")

Important_Features <- read.csv("Interpretation_Features_ED_LD.csv")
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_light(base_size = 16) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for interpretation") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 16))
plot_

#Variable selection procedure for prediction.
vsurf.pred <- VSURF_pred(Outcome~., data=mydata2, err.interp = vsurf.interp$err.interp, varselect.interp = vsurf.interp$varselect.interp)
plot(vsurf.pred, var.names = TRUE,
     nvar.pred = length(vsurf.pred$varselect.pred))

#use ggplot2 to plot our set of important features for prediction:
Important_Features <- data.frame(Feature = vsurf.pred$varselect.pred, Importance = vsurf.pred$err.pred)
write.csv(Important_Features, "Prediction_Features_ED_LD.csv")

ortho <- read.csv("Prediction_Features_ED_LD.csv")
cluster <- read.csv("Book1.csv")
idx <- match(ortho$Feature, cluster$Feature )
ortho$Biomarker <- cluster$Biomarker[ idx ]
write.csv(ortho, "Prediction_Features_ED_LD.csv")

Important_Features <- read.csv("Prediction_Features_ED_LD.csv")
plot_ <- ggplot(Important_Features, 
                aes(x= reorder(Feature,
                               (Importance)) , y = Importance, fill = Category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw(base_size = 12) +
  xlab("") + 
  ylab("Mean Decrease Accuracy (OOB error)") +
  ggtitle("Variable selection for prediction") + scale_fill_manual("legend", values = c("Clinical" = "black", "PB_biomarker" = "#800080")) +
  theme(axis.text.y = element_text(size = 12))
plot_

############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 5- Run random forest rmodel with the variables important for prediction as determined above - only data from days 0 and 3 and including clinical and luminex panel 1 data from patients d0

# 5.1- Recovered vs fatal

# train model
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_d0_d3_vsurf_v2.csv")
mydata <- multiplex[,-1]
mydata2 <- mydata[,1:85]
mydata2 <- as.data.frame(mydata2)
mydata2$Outcome <- as.factor(mydata2$Outcome)

#missing data in RF: https://stackoverflow.com/questions/8370455/how-to-use-random-forests-in-r-with-missing-values
set.seed(222)
mydata2.imputed <- rfImpute(Outcome ~ ., mydata2, iter=10, ntree=1000)
write.csv(mydata2.imputed, "COVID_Clinical_Luminex_deceased_recovered_d0_d3_vsurf_v2_imputed.csv")
mydata2.imputed <- read.csv("COVID_Clinical_Luminex_deceased_recovered_d0_d3_vsurf_v2_imputed.csv")
mydata2.imputed$Outcome <- as.factor(mydata2.imputed$Outcome)

set.seed(123)
sample_data <- initial_split(mydata2.imputed, prop=0.7, strata = "Outcome")
train_data <- training(sample_data)
test_data <- testing(sample_data)


set.seed(1) #running RF with imputed missing data increased massively the OOB error rate, so just try it for the predictions of disease progression
rf_train <- randomForest(Outcome~IL12p70+CCL1+IFNlambda2+PD_L1+IL33+IFNb+Creatinine+Platelets, data=train_data, importance=T, ntree=1000, nodesize=5)   #it does not handle missing data
print(rf_train)
plot(rf_train)
reprtree:::plot.getTree(rf_train, k=7, depth=4)# summarize the fitsummary(rf_train)

# make predictions in the test data set
predictions <- predict(rf_train, test_data)
# summarize accuracy and predicitons in the test data set
table(predictions, test_data$Outcome)
mean(predictions == test_data$Outcome)
confusionMatrix(predictions, test_data$Outcome)

test_plot<-ggplot(train_data, aes(IL12p70, IFNb, colour=Category)) + geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", y="IFNbeta (pg/mL)", x="IL12p70 (pg/mL)") + 
  theme_bw(base_size = 24) +  
  scale_color_manual(labels = c("Recovered", "Fatal"), values=c("dark blue", "light grey")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 4.8, linetype="dotted", linewidth = 1) +
  geom_vline(xintercept = 50, linetype="dotted", linewidth = 1)
test_plot

test_plot<-ggplot(train_data, aes(PD_L1, Creatinine, colour=Category)) +geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", y="Creatinine (mg/dL)", x="PD-L1 (pg/mL)") + 
  theme_bw(base_size = 24) + 
  scale_color_manual(labels = c("Fatal", "Recovered"), values=c("dark blue", "light grey")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 1.5, linetype="dotted", linewidth = 1) +
  geom_vline(xintercept = 270, linetype="dotted", linewidth = 1)
test_plot

test_plot<-ggplot(train_data, aes(PD_L1, IFNb, colour=Category)) + geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", y="IFNbeta") + 
  theme_bw(base_size = 24) +  
  scale_color_manual(values=c("dark blue", "light grey")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 4.8, linetype="dotted", linewidth = 1) +
  geom_vline(xintercept = 280, linetype="dotted", linewidth = 1)
test_plot

#performance of the random forests model applied in the train dataset
library(ROCR)
pred=predict(rf_train, type = "prob")
pred1 = prediction(pred[,2], train_data$Outcome)
# 1. True Positive and Negative Rate
pred2 = performance(pred1, "tpr", "fpr")
# 2. Plot the ROC curve
plot(pred2,main="ROC Curve for Random Forest",col="black",lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="black")
# 3. Area under curve
auc <- performance(pred1, measure = "auc")
auc
auc@y.values
# 4 - Accuracy
acc.perf = performance(pred1, measure = "acc")
plot(acc.perf)

#plotting with pROC and ggroc
roc1 <- multiclass.roc(train_data$Outcome, pred)
rs_Recov_Fatal <- roc1[['rocs']]
plot.roc(rs_Recov_Fatal[["0/1"]][[1]], legacy.axes=TRUE)

Recov_vs_Fatal <- rs_Recov_Fatal[["0/1"]][[1]]

ggroc(list(Recov_vs_Fatal_ROC=Recov_vs_Fatal), alpha = 1, size = 2, legacy.axes=TRUE) + theme_bw(base_size = 14) + ggtitle("ROC curve for RF") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype="dashed") + xlab("1-Specificity") + ylab("Sensitivity") +
  scale_colour_manual(values = c("dark blue")) + theme(legend.position="bottom", legend.text = element_text(size = 10))


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

#performance of the random forests model applied in the test dataset
library(ROCR)
pred1=predict(rf_train, newdata = test_data,type = "prob")
pred = prediction(pred1[,2], test_data$Outcome)
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

############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 5.2- Recovered vs Early death vs Late death - only data from days 0 and 3

# train model
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_d0_d3_vsurf_v2.csv")
mydata <- multiplex[,-1]
mydata2 <- mydata[,1:85]
mydata2 <- as.data.frame(mydata2)
mydata2$Outcome <- as.factor(mydata2$Outcome)

#missing data in RF: https://stackoverflow.com/questions/8370455/how-to-use-random-forests-in-r-with-missing-values
set.seed(222)
mydata2.imputed <- rfImpute(Outcome ~ ., mydata2, iter=10, ntree=1000)
write.csv(mydata2.imputed, "COVID_Clinical_Luminex_deceased_recovered_d0_d3_vsurf_v2_imputed.csv")

set.seed(123)
sample_data <- initial_split(mydata2.imputed, prop=0.7, strata = "Outcome")
train_data <- training(sample_data)
test_data <- testing(sample_data)

set.seed(1)
rf_train <- randomForest(Outcome~IL12p70+IFNlambda2+IFNb+IL33+CCL1+PD_L1+NCLR+ALT, data=train_data, importance=T, ntree=1000, nodesize=5)   #mtry=5, ntree=1000
print(rf_train)
plot(rf_train)
reprtree:::plot.getTree(rf_train, k=2, depth = 4)
# summarize the fit
summary(rf_train)

# make predictions in the test data set
predictions <- predict(rf_train, test_data)
# summarize accuracy and predicitons in the test data set
table(predictions, test_data$Outcome)
mean(predictions == test_data$Outcome)
confusionMatrix(predictions, test_data$Outcome)

test_plot<-ggplot(train_data, aes(IL12p70, PD_L1, colour=Outcome)) + geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", y="PD-L1") + 
  theme_light(base_size = 16) +  
  scale_color_manual(values=c("dark grey", "salmon","light blue")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 600, linetype="dotted") +
  geom_vline(xintercept = 30, linetype="dotted")
test_plot

test_plot<-ggplot(train_data, aes(ALT, PD_L1, colour=Outcome)) + geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", y="PD-L1") + 
  theme_light(base_size = 16) +  
  scale_color_manual(values=c("dark grey", "salmon","light blue")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 600, linetype="dotted") +
  geom_vline(xintercept = 80, linetype="dotted")
test_plot

#performance of the random forests model applied 
library(pROC) # multiclass: Multi-class AUC
library(ROCR)
# AUROC in the train dataset
pred=predict(rf_train, type = "prob")
roc1 <- multiclass.roc(train_data$Outcome, pred, levels = c(0, 1)) # AUROC in the train dataset
roc1
# AUROC in the test dataset
pred1=predict(rf_train, newdata=test_data,type = "prob")
roc2 <- multiclass.roc(test_data$Outcome, pred1)
roc2

rs <- roc1[['rocs']]
plot.roc(rs[["1/2"]][[1]])

Recov_vs_ED <- rs[["0/1"]][[1]]
Recov_vs_LD <- rs[["0/2"]][[1]]
ED_vs_LD <- rs[["1/2"]][[1]]

ggroc(list(ED_vs_LD_ROC=ED_vs_LD, Recov_vs_ED_ROC=Recov_vs_ED, Recov_vs_LD_ROC= Recov_vs_LD), alpha = 1, size = 2) + theme_bw(base_size = 14) + ggtitle("ROC curve for RF") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed") + xlab("False positive rate") + ylab("True positive rate") +
  scale_colour_manual(values = c("salmon", "light blue", "grey")) + theme(legend.position="bottom", legend.text = element_text(size = 10))

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

#####################################################################################################################################################################################################################################################################################################################################################################

# 5.3- Early death vs Late death - only data from days 0 and 3

# train model
multiplex <- read.csv("./RF_data_d0_d3_missing_data/COVID_Clinical_Luminex_ED_LD_d0_d3_vsurf_v2.csv")
mydata <- multiplex[,-1]
mydata2 <- mydata[,1:85]
mydata2 <- as.data.frame(mydata2)
mydata2$Outcome <- as.factor(mydata2$Outcome)

#missing data in RF: https://stackoverflow.com/questions/8370455/how-to-use-random-forests-in-r-with-missing-values
set.seed(222)
mydata2.imputed <- rfImpute(Outcome ~ ., mydata2, iter=10, ntree=1000)
write.csv(mydata2.imputed, "COVID_Clinical_Luminex_ED_LD_d0_d3_vsurf_v2_imputed.csv")
mydata2.imputed <- read.csv("COVID_Clinical_Luminex_ED_LD_d0_d3_vsurf_v2_imputed.csv")

set.seed(123)
sample_data <- initial_split(mydata2.imputed, prop=0.7, strata = "Outcome")
train_data <- training(sample_data)
test_data <- testing(sample_data)

# RF model if using VSURF without imputed missing data
# set.seed(1234)
# rf_train <- randomForest(Outcome~PD_L1+CCL1, data=train_data, importance=T, ntree=1000)   #mtry=5, ntree=1000
# print(rf_train)
# plot(rf_train)
# reprtree:::plot.getTree(rf_train, k=70, depth = 4)
# summarize the fit
# summary(rf_train)


# RF model if using VSURF from imputed missing data - this is the one for the paper!
set.seed(1234)
rf_train <- randomForest(Outcome~PD_L1+CCL1+LDH+Creatinine, data=train_data, importance=T, ntree=1000)   #mtry=5, ntree=1000
print(rf_train)
plot(rf_train)
reprtree:::plot.getTree(rf_train, k=800, depth = 4)# summarize the fitsummary(rf_train)
# make predictions in the test data set
predictions <- predict(rf_train, test_data)
# summarize accuracy and predicitons in the test data set
table(predictions, test_data$Outcome)
mean(predictions == test_data$Outcome)
confusionMatrix(predictions, test_data$Outcome)

test_plot<-ggplot(train_data, aes(Urea, LDH, colour=Category2)) + geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set",  x="Creatinine (mg/dL)", y="LDH (U/L)") + 
  theme_bw(base_size = 24) +  
  scale_color_manual(labels = c("Early death" ,"Late death"), values=c("salmon", "light blue")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 1200, linetype="dotted", linewidth = 1) +
  geom_vline(xintercept = 10, linetype="dotted", linewidth = 1)
test_plot

test_plot<-ggplot(train_data, aes(PD_L1, CCL1, colour=Category2)) + geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", x="PD_L1 (pg/mL)", y="CCL1 (pg/mL)") + 
  theme_bw(base_size = 24) +  
  scale_color_manual(labels = c("Early death" ,"Late death"), values=c("salmon", "light blue")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 20, linetype="dotted", linewidth = 1) +
  geom_vline(xintercept = 600, linetype="dotted", linewidth = 1)
test_plot

test_plot<-ggplot(train_data, aes(M_CSF, CXCL2, colour=Outcome)) + geom_jitter(size=8) +
  labs(title="Distribution of Fatal COVID-19 Outcome in Test Set", x="M_CSF (pg/mL)", y="CXCL2 (pg/mL)") + 
  theme_bw(base_size = 16) +  
  scale_color_manual(labels = c("Late death", "Early death"), values=c("light blue", "salmon")) +  
  theme(legend.position="bottom") +
  geom_hline(yintercept = 400, linetype="dotted") +
  geom_vline(xintercept = 3000, linetype="dotted")
test_plot

#performance of the random forests model applied in the train dataset
library(ROCR)
pred=predict(rf_train, type = "prob")
pred1 = prediction(pred[,2], train_data$Outcome)
# 1. True Positive and Negative Rate
pred2 = performance(pred1, "tpr","fpr")
# 2. Plot the ROC curve
plot(pred2,main="ROC Curve for Random Forest",col="black",lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="black")
# 3. Area under curve
auc <- performance(pred1, measure = "auc")
auc
auc@y.values
# 4 - Accuracy
acc.perf = performance(pred1, measure = "acc")
plot(acc.perf)

#plotting with pROC and ggroc
roc1 <- multiclass.roc(train_data$Outcome, pred)
rs <- roc1[['rocs']]
plot.roc(rs[["0/1"]][[1]])

ED_vs_LD <- rs[["0/1"]][[1]]

ggroc(list(ED_vs_LD_ROC=ED_vs_LD), alpha = 1, size = 2) + theme_bw(base_size = 14) + ggtitle("ROC curve for RF") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed") + xlab("False positive rate") + ylab("True positive rate") +
  scale_colour_manual(values = c("dark blue")) + theme(legend.position="bottom", legend.text = element_text(size = 10))

#plotting Recovery vs Fatal and ED vs LD curves with pROC and ggroc
ggroc(list(Recov_vs_Fatal_ROC=Recov_vs_Fatal, ED_vs_LD_ROC=ED_vs_LD), alpha = 1, size = 2) + theme_bw(base_size = 14) + ggtitle("ROC curve for RF") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="dashed") + xlab("False positive rate") + ylab("True positive rate") +
  scale_colour_manual(values = c("dark blue", "red")) + theme(legend.position="bottom", legend.text = element_text(size = 10))


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

#performance of the random forests model applied in the test dataset
library(ROCR)
pred1=predict(rf_train, newdata = test_data,type = "prob")
pred = prediction(pred1[,2], test_data$Outcome)
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

