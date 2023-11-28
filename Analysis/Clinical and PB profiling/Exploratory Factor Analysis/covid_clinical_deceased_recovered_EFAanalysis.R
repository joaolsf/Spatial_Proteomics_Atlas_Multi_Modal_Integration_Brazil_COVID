library(psych)
library(pheatmap)
library(ggplot2)
library(corrplot)
library(car)


#1- Data Pre-processing
# reading file (importar o arquivo)
multiplex <- read.csv("COVID_Clinical_deceased_recovered_EFA.csv")

# remove columns
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c( "Category", "Category2", "Day_of_Sampling", "Average_Day_of_Symptoms", "Day_of_Symptoms","DOS2")]
multiplex <- multiplex[,select_cols]

df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

# Prepare data
mydata <- scale(df_to_cluster)

#2- The Factorability of the Data
X <- data[,-c(13)]
Y <- data[,13]

#2.1- KMO
#The Kaiser-Meyer-Olkin (KMO) used to measure sampling adequacy is a better measure of factorability.
KMO(r=cor(mydata))
#According to Kaiser’s (1974) guidelines, a suggested cutoff for determining the factorability of the sample data is KMO ≥ 60.
#The total KMO is 0.65, indicating that, based on this test, we can probably conduct a factor analysis.

#2.2- Bartlett’s Test of Sphericity
cortest.bartlett(mydata)
#Small values (8.84e-290 < 0.05) of the significance level indicate that a factor analysis may be useful with our data.
det(cor(mydata))
#We have a positive determinant, which means the factor analysis will probably run.

#3- The Number of Factors to Extract
#3.1- Scree Pilot
fafitfree <- fa(mydata,nfactors = ncol(mydata), rotate = "none")

n_factors <- length(fafitfree$e.values)

scree <- data.frame(Factor_n =  as.factor(1:n_factors), Eigenvalue = fafitfree$e.values)

ggplot(scree, aes(x = Factor_n, y = Eigenvalue, group = 1)) + 
  geom_point() + geom_line() +
  xlab("Number of factors") +
  ylab("Initial eigenvalue") +
  labs(title = "Scree Plot", 
        subtitle = "(Based on the unreduced correlation matrix)")

#3.2- Parallel Analysis
#We can use the parallel() function from the nFactors package (Raiche & Magis, 2020) to perform a parallel analysis.
parallel <- fa.parallel(mydata)

#4- Conducting the Factor Analysis
#4.1- Factor analysis using fa method
fa.none <- fa(r=mydata, 
              nfactors = 3, 
              # covar = FALSE, SMC = TRUE,
              fm="pa", # type of factor analysis we want to use (“pa” is principal axis factoring)
              max.iter=100, # (50 is the default, but we have changed it to 100
              rotate="varimax") # none rotation
print(fa.none)

#4.2- Graph Factor Loading Matrices
fa.diagram(fa.none, simple=TRUE, side=2, cex=1,marg=c(.2,8,.2,8),adj=4,ic=FALSE, cut=.1)
fa.graph(fa.none)
head(fa.none$scores)
head(fa.none$loadings)

write.csv(fa.none$scores, "COVID_Clinical_deceased_recovered_factor_scores_iteration2.csv")
write.csv(fa.none$loadings, "COVID_Clinical_deceased_recovered_factor_loadings_iteration2.csv")

fa.none
regdata <- cbind(mydata["QD"], fa.none$scores)
#Labeling the data
names(regdata) <- c("F1", "F2",
                     "F3", "F4", "F5")
head(regdata)

#4.3- Factor analysis using the factanal method
factanal.none <- factanal(mydata, factors=6, scores = c("regression"), rotation = "varimax")
print(factanal.none)




