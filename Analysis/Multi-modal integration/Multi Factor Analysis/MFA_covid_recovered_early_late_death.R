#Multiple Factor Analysis 

library("FactoMineR")
library("factoextra")
library("missMDA")
library(RColorBrewer)
library(ggrepel)
library(viridis)
library(ggpubr)
options(ggrepel.max.overlaps = 100)

multiplex <- read.csv("COVID_Clinical_Luminex_IMC_recovered_early_late_death_MFA.csv")
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("BlockID", "Length","Score", "Age", "Weight", "BMI", "DOSAA", "DDOUO","DOSAPA", "DOS_original")]
multiplex <- multiplex[,select_cols]
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

#subsetting active variables for the PCA
#df_to_cluster.active <- df_to_cluster[,1:232]
#df_to_cluster.inactive <- df_to_cluster[,233:235]

# Prepare input data
mydata <- (df_to_cluster)

#remove Recovered patients
mydata <- mydata[1:74,]

mydata_luminex <- mydata[,1:72]
mydata_clinical <- mydata[,73:88]
mydata_IMC <- mydata[,89:201]
mydata_sup <- mydata[,202:204]

#1-Handling missing data in the dataframe:
#select the number of dimensions that will be used in the algorithm using the function estim_ncpPCA
#Perform the (regularized) iterative MFA algorithm with the number of dimensions selected in the previous step, using the function imputePCA:
mydata2 <- imputeMFA(mydata[,1:201], group = c(72, 16, 113), #luminex, clinical, IMC+histo continuous, histo categorical
                     type = c("s", "s", "s"),
                     ncp = 5)

#perform MFA on the imputed data set. To this end, we propose to use the MFAfunction of the FactoMineR package:
#The function MFA()[FactoMiner package] can be used. A simplified format is :
#MFA (base, group, type = rep("s",length(group)), ind.sup = NULL, name.group = NULL, num.group.sup = NULL, graph = TRUE)

#The R code below performs the MFA on the wines data using the groups: odor, visual, odor after shaking and taste. 
#These groups are named active groups. The remaining group of variables - origin (the first group) and overall judgement (the sixth group) - are named supplementary groups; 
#num.group.sup = c(1, 6):
mydata2 <- cbind.data.frame(mydata2$completeObs, mydata_sup)

#write.csv(mydata2, "COVID_Clinical_Luminex_IMC_recovered_early_late_death_MFA_Complete.csv")
#mydata <- read.csv("Complete.csv")
# data frame para cluster
#mydata2 <- mydata[,-1]
#rownames(mydata2) <- mydata$RecordID

mydata2[,202:204] <- lapply(mydata2[,202:204],factor)

mydata3 <- MFA(mydata2, group = c(72, 16, 113, 1, 1, 1), #luminex, clinical, IMC+histo continuous, histo categorical
               type = c("s", "s", "s", "n", "n", "n"),
               name.group = c("PB_Luminex","Clinical","IMC",
                              "Category", "Day_of_sampling", "Days_of_symptoms"), 
               num.group.sup = c(4, 5, 6),
               graph = FALSE)
print(mydata3)

#Visualization and interpretation
#We’ll use the factoextra R package to help in the interpretation and the visualization of the multiple factor analysis.
#The functions below [in factoextra package] will be used:
#get_eigenvalue(res.mfa): Extract the eigenvalues/variances retained by each dimension (axis).
#fviz_eig(res.mfa): Visualize the eigenvalues/variances.
#get_mfa_ind(res.mfa): Extract the results for individuals.
#get_mfa_var(res.mfa): Extract the results for quantitative and qualitative variables, as well as, for groups of variables.
#fviz_mfa_ind(res.mfa), fviz_mfa_var(res.mfa): Visualize the results for individuals and variables, respectively.

#1- Eigenvalues / Variances
#The proportion of variances retained by the different dimensions (axes) can be extracted using the function get_eigenvalue() [factoextra package] as follow:
eig.val <- get_eigenvalue(mydata3)
head(eig.val)
#The function fviz_eig() or fviz_screeplot() [factoextra package] can be used to draw the scree plot:
fviz_screeplot(mydata3)

#2- Graph of variables
#2.1- Groups of variables
#The function get_mfa_var() [in factoextra] is used to extract the results for groups of variables.
#This function returns a list containing the coordinates, the cos2 and the contribution of groups.
group <- get_mfa_var(mydata3, "group")
group

#The different components can be accessed as follow: (check the PCA R files with this package to see how to plot the parameters below using corrplot)
# Coordinates of groups
head(group$coord)
# Cos2: quality of representation on the factore map
head(group$cos2)
# Contributions to the  dimensions
head(group$contrib)

#To plot the groups of variables, type this:
#red color = active groups of variables
#green color = supplementary groups of variables
fviz_mfa_var(mydata3, "group")

#To draw a bar plot of groups contribution to the dimensions, use the function fviz_contrib():
# Contribution to the first dimension (check the PCA R files with this package to see how to edit this plot to include more parameters)
fviz_contrib(mydata3, "group", axes = 1)
# Contribution to the second dimension
fviz_contrib(mydata3, "group", axes = 2)

#2.2- Quantitative variables
#The function get_mfa_var() [in factoextra] is used to extract the results for quantitative variables. 
#This function returns a list containing the coordinates, the cos2 and the contribution of variables:
quanti.var <- get_mfa_var(mydata3, "quanti.var")
quanti.var 

#The different components can be accessed as follow:
# Coordinates
head(quanti.var$coord)
# Cos2: quality on the factor map
head(quanti.var$cos2)
# Contributions to the dimensions
head(quanti.var$contrib)

#In this section, we’ll describe how to visualize quantitative variables colored by groups.
#Next, we’ll highlight variables according to either:
#i) their quality of representation on the factor map or 
#ii) their contributions to the dimensions.

#Correlation between quantitative variables and dimensions. The R code below plots quantitative variables colored by groups. 
#The argument palette is used to change group colors (see ?ggpubr::ggpar for more information about palette). 
#Supplementary quantitative variables are in dashed arrow and violet color. We use repel = TRUE, to avoid text overlapping.
fviz_mfa_var(mydata3, "quanti.var", 
             palette = "jco", 
             repel = TRUE, 
             select.var = list(contrib = 60))

#To make the plot more readable, we can use geom = c(“point”, “text”) instead of geom = c(“arrow”, “text”). 
#We’ll change also the legend position from “right” to “bottom”, using the argument legend = “bottom”:
fviz_mfa_var(mydata3, 
             "quanti.var", 
             palette = c("#F79C79", "#21908CFF", "#440154FF"), 
             repel = TRUE, 
             select.var = list(cos2 = 60),
             geom = c("point", "text"),
             shape.var=19,
             legend = "bottom")

fviz_mfa_var(mydata3, 
             "quanti.var", 
             palette = c("#F79C79", "#21908CFF", "#440154FF"), 
             repel = TRUE,
             select.var = list(contrib = 140),
             geom = c("point", "text"),
             shape.var= 19, labelsize = 2,
             legend = "bottom")  + theme(text = element_text(size = 5),
                                      axis.title = element_text(size = 5),
                                      axis.text = element_text(size = 5)
             )

#The contribution of quantitative variables (in %) to the definition of the dimensions can be visualized using the function fviz_contrib() [factoextra package]. Variables are colored by groups. The R code below shows the top 20 variable categories contributing to the dimensions:
# Contributions to dimension 1
fviz_contrib(mydata3, choice = "quanti.var", axes = 1, top = 60, xtickslab.rt = 90,
             palette = c("#F79C79", "#21908CFF", "#440154FF"), ggtheme = theme_light())
# Contributions to dimension 2
fviz_contrib(mydata3, choice = "quanti.var", axes = 2, top = 60,xtickslab.rt = 90,
             palette = c("#F79C79", "#21908CFF", "#440154FF"), ggtheme = theme_light())

#The most contributing quantitative variables can be highlighted on the scatter plot using the argument col.var = “contrib”. 
#This produces a gradient colors, which can be customized using the argument gradient.cols.
fviz_mfa_var(mydata3, "quanti.var", col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE)
#Similarly, you can highlight quantitative variables using their cos2 values representing the quality of representation on the factor map. 
#If a variable is well represented by two dimensions, the sum of the cos2 is closed to one. 
#For some of the row items, more than 2 dimensions might be required to perfectly represent the data.
# Color by cos2 values: quality on the factor map
fviz_mfa_var(mydata3, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             col.var.sup = "violet", repel = TRUE)
#To create a bar plot of variables cos2, type this:
fviz_cos2(mydata3, choice = "quanti.var", axes = 1)

#3- Graph of individuals
#To get the results for individuals, type this:
ind <- get_mfa_ind(mydata3)
ind
#To plot individuals, use the function fviz_mfa_ind() [in factoextra].
#By default, individuals are colored in blue. 
#However, like variables, it’s also possible to color individuals by their cos2 values:
fviz_mfa_ind(mydata3, col.ind = "contrib", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE)
#In the plot above, the supplementary qualitative variable categories are shown in black. Env1, Env2, Env3 are the categories of the soil. Saumur, Bourgueuil and Chinon are the categories of the wine Label. 
#If you don’t want to show them on the plot, use the argument invisible = “quali.var”

#Note that, it’s possible to color the individuals using any of the qualitative variables in the initial data table. 
#To do this, the argument habillage is used in the fviz_mfa_ind() function. 
#For example, if you want to color the wines according to the supplementary qualitative variable “Label”, type this:
p <- fviz_mfa_ind(mydata3,
             geom.ind = c("point"),
             repel = TRUE, 
             shape.ind = 19,
             pointsize = 5,
             habillage = "Day",
             fill.ind = "Day",
             palette = magma(6, alpha = 1, begin = 0, end = 1, direction = -1),
             col.ind = "white", 
             label = "var",
             addEllipses = FALSE, ellipse.type = "confidence", invisible="quali") + theme_bw()

p + geom_point(aes(shape = mydata$Category), size=5) + scale_shape_manual(values=c(21, 22))

fviz_mfa_ind(mydata3,
             geom.ind = c("point"),
             repel = TRUE, 
             shape.ind = 19,
             pointsize = 5,
             habillage = "Category",
             fill.ind = "Category",
             palette = c("salmon", "light blue"),
             col.ind = "white", 
             label = "var",
             addEllipses = TRUE, ellipse.type = "confidence", invisible="quali") + theme_bw()


#If you want to color individuals using multiple categorical variables at the same time, use the function fviz_ellipses() [in factoextra] as follow:
fviz_ellipses(mydata3, c("Label", "Soil"), repel = TRUE)
#Alternatively, you can specify categorical variable indices:
fviz_ellipses(mydata3, 1:2, geom = "point")

#3.1- Graph of partial individuals
#The results for individuals obtained from the analysis performed with a single group are named partial individuals. 
#In other words, an individual considered from the point of view of a single group is called partial individual.

#In the default fviz_mfa_ind() plot, for a given individual, the point corresponds to the mean individual or the center of gravity of the partial points of the individual. 
#That is, the individual viewed by all groups of variables.
#For a given individual, there are as many partial points as groups of variables.
#The graph of partial individuals represents each wine viewed by each group and its barycenter. 
#To plot the partial points of all individuals, type this:
fviz_mfa_ind(mydata3, partial = "all") 
#If you want to visualize partial points for wines of interest, let say c(“1DAM”, “1VAU”, “2ING”), use this:
fviz_mfa_ind(mydata3a, partial = c("1DAM", "1VAU", "2ING")) 

#4- Graph of partial axes
#The graph of partial axes shows the relationship between the principal axes of the MFA and the ones obtained from analyzing each group using either a PCA (for groups of continuous variables) or a MCA (for qualitative variables).
fviz_mfa_axes(mydata3)

#It can be seen that, he first dimension of each group is highly correlated to the MFA’s first one. The second dimension of the MFA is essentially correlated to the second dimension of the olfactory groups.
