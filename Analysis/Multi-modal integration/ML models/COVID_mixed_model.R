
# Save an object to a file
saveRDS(covid_clinical, file = "covid_clinical.rds")
# Restore the object
readRDS(file = "my_data.rds")

# Nested Mixed models
# We want to use all the data, but account for the data coming from Age and gender
# let's add Age and gender as a fixed effect to our model

library(lme4)

covid.lm <- lmer(Lymphocytes ~ Age + Gender + (1|patient), data = covid_clinical)

summary(covid.lm)

library(stargazer)
stargazer(covid.lm, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")


#As always, it’s good practice to have a look at the plots to check our assumptions:
plot(covid.lm)
qqnorm(resid(covid.lm))
qqline(resid(covid.lm))

library(ggplot2)
(mm_plot <- ggplot(covid_clinical, aes(x = Age, y = Hemmoglobin, colour = patient)) +
    facet_wrap(~day, nrow=2) +   # a panel for each mountain range
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(covid_clinical, pred = predict(covid.lm)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)


library(ggeffects)  # install the package first if you haven't already, then load it

# Extract the prediction data frame
pred.mm1 <- ggpredict(covid.lm, terms = c("Age"))  # this gives overall predictions for the model
pred.mm2 <- ggpredict(covid.lm, terms = c("Gender"))

# Plot the predictions 
(ggplot(pred.mm1) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = covid_clinical,                      # adding the raw data (scaled values)
               aes(x = Age, y = Lymphocytes, colour = day)) + 
    labs(x = "Age (indexed)", y = "Lymphocytes", 
         title = "Age does not affect lymphocytes in patients") + 
    theme_minimal()
)

(ggplot(pred.mm2) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = covid_clinical,                      # adding the raw data (scaled values)
               aes(x = Gender, y = Urea, colour = day)) + 
    labs(x = "Sex", y = "Urea", 
         title = "Sex does not affect urea in patients") + 
    theme_minimal()
)

(covid.p <- ggplot(covid_clinical, aes(x = Age, y = Hemoglobin)) +
    geom_point(aes(colour = Gender)) +                                # scatter plot, coloured by sex
    labs(x = "Days of Hospital Admission", y = "Hemoglobin") +
    stat_smooth(method = "lm", aes(fill = Gender, colour = Gender)) +    # adding regression lines for each sex
    scale_colour_manual(values = c("#FFC125", "#36648B")) +
    scale_fill_manual(values = c("#FFC125", "#36648B")) +
    theme_minimal())

# Other useful functions 
coefficients(covid.lm) # model coefficients
confint(covid.lm, level=0.95) # CIs for model parameters 
predicted <- fitted(covid.lm) # predicted values - original value minus residuals
res <- residuals(covid.lm) # residuals
anova(covid.lm) # anova table 
vcov(covid.lm) # covariance matrix for model parameters 
influence(covid.lm) # regression diagnostics

write.csv((broom::tidy(predicted)), "Urea_sex_adj_values.csv")
write.csv((broom::tidy(res)), "Urea_sex_residuals.csv")

# K-fold cross-validation
library(DAAG)
cv.lm(df=covid.lm, fit, m=10) # 10 fold cross-validation
write.csv((broom::tidy(lm_results)), "lm_results_age_gender.csv")


library(ggplot2)
ggplot(covid_clinical, aes(x = Gender, y = Hemoglobin)) +
  geom_point()+
  geom_smooth(method = "lm")


lm_results <- lm(Respuratory.rate ~ Age + Gender, data = covid_clinical)
summary(lm_results)

predicted <- fitted(lm_results) # predicted values - original value minus residuals
res <- residuals(lm_results) # residuals
write.csv((broom::tidy(predicted)), "hemoglobin_gender_adj_values.csv")
write.csv((broom::tidy(res)), "hemoglobin_gender_residuals.csv")

# match covid_data with adjusted values data frame by names
covid_data_adjusted <- match(covid_clinical, Age_gender_adjusted_values,by=c("names", "names.2", "names.3", "names.4"))

### Assumptions?

## Plot the residuals - the red line should be close to being flat, like the dashed grey line

plot(covid.lm, which = 1)  # not perfect, but look alright

## Have a quick look at the  qqplot too - point should ideally fall onto the diagonal dashed line

plot(covid.lm, which = 2)  # a bit off at the extremes, but that's often the case; again doesn't look too bad

## Have a look at the data to see if above is true
boxplot(y ~ Age, data = clinical_data)  # certainly looks like something is going on here

## We could also plot it colouring points by gender
ggplot(clinical_data, aes(x = Age, y = Hemoglobin, colour = Gender))+
  geom_point(size = 2)+
  theme_classic()+
  theme(legend.position = "none")

## From the above plots it looks like our mountain ranges vary both in the dragon body length and in their test scores. This confirms that our observations from within each of the ranges aren't independent. We can't ignore that.

## So what do we do?

###----- Run multiple analyses -----###

## We could run many separate analyses and fit a regression for each of the mountain ranges.

## Lets have a quick look at the data split by mountain range
## We use the facet_wrap to do that

ggplot(aes(bodyLength, testScore), data = dragons) + geom_point() +
  facet_wrap(~ mountainRange) +
  xlab("length") + ylab("test score")

#merge for downstream analysis

ortho <- read.csv("covid_data.csv")
cluster1 <- read.csv("Hemoglobin_adjusted.csv")

idx <- match(ortho$names, cluster1$names )
ortho$Hemoglobin_residuals <- cluster1$Hemoglobin_residuals [ idx ]
ortho$Hemoglobin_predicted <- cluster1$Hemoglobin_predicted [ idx ]

cluster2 <- read.csv("Lymphocyte_adjusted.csv")

idx2 <- match(ortho$names, cluster2$names )
ortho$Lymphocyte_residuals <- cluster2$Lymphocyte_residuals [ idx2 ]
ortho$Lymphocyte_predicted <- cluster2$Lymphocyte_predicted [ idx2 ]

cluster3 <- read.csv("Platelet_adjusted.csv")

idx3 <- match(ortho$names, cluster3$names )
ortho$Platelet_residuals <- cluster3$Platelet_residuals [ idx3 ]
ortho$Platelet_predicted <- cluster3$Platelets_predicted [ idx3 ]

cluster4 <- read.csv("Urea_adjusted.csv")

idx4 <- match(ortho$names, cluster4$names )
ortho$Urea_residuals <- cluster4$Urea_residuals [ idx4 ]
ortho$Urea_predicted <- cluster4$Urea_predicted [ idx4 ]

write.csv(ortho, file="Age_sex_adjusted_values.csv")

df <- read.csv("Age_sex_adjusted_values.csv")
library(tidyr)
df2 <- df[order(df$Annotation),]
write.csv(df2, file="Age_sex_adjusted_values_organized.csv")

ortho <- read.csv("RecordID.csv")
cluster5 <- read.csv("Age_sex_adjusted_values.csv")
idx4 <- match(ortho$day, cluster5$day)
idx5 <- match(ortho$Record.ID, cluster5$PatientID)

ortho$Hemoglobin_residuals <- cluster5$Hemoglobin_residuals [ idx5 ]
ortho$Hemoglobin_predicted <- cluster5$Hemoglobin_predicted [ idx5 ]
ortho$Lymphocyte_residuals <- cluster5$Lymphocyte_residuals [ idx5 ]
ortho$Lymphocyte_predicted <- cluster5$Lymphocyte_predicted [ idx5 ]
ortho$Platelet_residuals <- cluster5$Platelet_residuals [ idx5 ]
ortho$Platelet_predicted <- cluster5$Platelets_predicted [ idx5 ]
ortho$Urea_residuals <- cluster5$Urea_residuals [ idx5 ]
ortho$Urea_predicted <- cluster5$Urea_predicted [ idx5 ]

write.csv(ortho, file="Age_sex_adjusted_values_organized.csv")
