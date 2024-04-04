
multiplex <- read.csv("/Users/joaoluizsfilho/Dropbox/Work Files/Matthia's Lab/Projects/Severe Anemia/Integration/cohort_anemia_lm.csv")

select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("Study_number","Group","CXCL12", "PB_parasitemia", 
                                               "PB_density", "BM_asexual_dens", "BM_gam_dens", "BM_parasitemia",
                                               "Weight_for_age_z_score", "Height_for_age_z_score", "Weight_for_length_z_score",
                                               "Treatment_anemia", "Anemia", "Treatment_Platelets", "Treatment_RD", "Treatment_CC",
                                               "Treatment_Tachycardia", "Treatment_chronic_wasting",
                                               "Treatment_chronic_stunting", "Luminex_cluster", "Clinical_cluster", "Clinical")]

# select_cols
multiplex <- multiplex[,select_cols]
# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

mydata <- df_to_cluster

model1 <- lm(mydata$Syndecan.1 ~ mydata$Age)
summary(model1)

model2 <- lm(mydata$Syndecan.1 ~ mydata$Age + mydata$Hemoglobin)
summary(model2)
confint(model2)
#The first step in interpreting the multiple regression analysis is to examine the F-statistic and the associated p-value, at the bottom of model summary.
#This means that, at least, one of the predictor variables is significantly related to the outcome variable.
#For a given the predictor, the t-statistic evaluates whether or not there is significant association between the predictor and the outcome variable, that is whether the beta coefficient of the predictor is significantly different from zero.
#For a given predictor variable, the coefficient (b) can be interpreted as the average effect on y of a one unit increase in predictor, holding all other predictors fixed.


Syndecan = 7351.176 - 9.231*mydata$Age - 479.914*mydata$Hemoglobin

#the overall quality of the model can be assessed by examining the R-squared (R2) and Residual Standard Error (RSE).
#R2 represents the proportion of variance, in the outcome variable y, that may be predicted by knowing the value of the x variables. An R2 value close to 1 indicates that the model explains a large portion of the variance in the outcome variable.
#A problem with the R2, is that, it will always increase when more variables are added to the model, even if those variables are only weakly associated with the response (James et al. 2014). A solution is to adjust the R2 by taking into account the number of predictor variables.
#The adjustment in the “Adjusted R Square” value in the summary output is a correction for the number of x variables included in the prediction model
#Residual Standard Error (RSE), or sigma: The RSE estimate gives a measure of error of prediction. The lower the RSE, the more accurate the model (on the data in hand).
#The error rate can be estimated by dividing the RSE by the mean outcome variable:
sigma(model1)/mean(mydata$Syndecan.1)

model3 <- lm(mydata$EPO ~ mydata$Age + mydata$Hemoglobin)
summary(model3)
sigma(model3)/mean(mydata$EPO)


y <- as.matrix(mydata[1:63])
lm_results <- lm(y ~ Age + Hemoglobin, data = mydata)
write.csv((broom::tidy(lm_results)), "lm_results_age.csv")

library(lm.beta)
regression_results <-broom::tidy(lm_results)
standardized_coefficients <- lm.beta(lm_results)
age_standardize_results <- coef(standardized_coefficients)

write.csv(age_standardize_results, "lm_results_age_standardized_coefficients.csv")

#Get all column names to run regression on
depVarList = setdiff(colnames(myData), c("date", "mktrf", "hml", "smb"))

#Loop over them and create model for each
allModels = lapply(depVarList, function(x){
  lm(formula= paste0("`", x, "` ~ mktrf + hml + smb"), 
     data= myData ,na.action = na.omit)
  
})

#Name the list of models to the column name
names(allModels) = depVarList

#Recuperate the columns of interest from the original dataframe
allResiduals = myData %>% select(-mktrf , -hml , -smb)

#Replace all values with residuals (ignore NA)
allResiduals[,-1] = sapply(2:ncol(allResiduals), function(x){
  residuals = allResiduals[,x]
  residuals[!is.na(residuals)] = allModels[[x-1]]$residuals
  residuals
})
