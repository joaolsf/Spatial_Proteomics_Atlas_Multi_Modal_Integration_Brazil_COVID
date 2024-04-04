
multiplex <- read.csv("COVID_Luminex_fatal_recovered_PCA.csv")
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

#subsetting active variables for the PCA
df_to_cluster.active <- df_to_cluster[,1:72]

# Prepare data
mydata <- scale(df_to_cluster.active)

#1- PCA Analysis using function PCA()[FactoMineR package]
library("FactoMineR")
library("factoextra")
library(RColorBrewer)
library(ggrepel)
library(viridis)
    
pca.mydata <- PCA(mydata, scale.unit = FALSE, ncp = 10, graph = TRUE)
print(pca.mydata)

#extract the eigenvalues/variances of PCs
eig.val <- get_eigenvalue(pca.mydata)
eig.val
#visualize the eigenvalues
fviz_eig(pca.mydata)
#extrcat the results for individuals and variables, respectively
get_pca_ind(pca.mydata)
get_pca_var(pca.mydata)
#Visualize the results individuals and variables, respectively.
fviz_pca_ind(pca.mydata)
fviz_pca_var(pca.mydata)
fviz_eig(pca.mydata, addlabels = TRUE, ylim = c(0, 50))

#method to extract the results, for variables,
var <- get_pca_var(pca.mydata)
var
fviz_pca_var(pca.mydata, col.var = "black")

#The quality of representation of the variables on factor map is called cos2 (square cosine, squared coordinates)
#visualize the cos2 of variables on all the dimensions using the corrplot package:
library("corrplot")
corrplot(var$cos2, method = "circle", is.corr=FALSE, tl.cex = 0.8, tl.col = "black", cl.ratio = 0.15, cl.align.text = "l")
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(pca.mydata, choice = "var", axes = 1:2)
# Color by cos2 values: quality on the factor map
fviz_pca_var(pca.mydata, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, alpha.var = "cos2") # Avoid text overlapping

#The contributions of variables in accounting for the variability in a given principal component are expressed in percentage.
#use the function corrplot() [corrplot package] to highlight the most contributing variables for each dimension:
corrplot(var$contrib, method = "circle", is.corr=FALSE, tl.cex = 0.8, tl.col = "black", cl.ratio = 0.15, cl.align.text = "l")
# Contributions of variables to PC1
fviz_contrib(pca.mydata, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(pca.mydata, choice = "var", axes = 2, top = 10)
#The total contribution to PC1 and PC2 is obtained with the following R code:
fviz_contrib(pca.mydata, choice = "var", axes = 1:5, top = 10)
# Color by contrib values: quality on the factor map
fviz_pca_var(pca.mydata, col.var = "contrib", 
             repel = TRUE, alpha.var = "contrib") + scale_color_gradient2(low="white", mid="blue", 
            high="red", midpoint=1) + theme_light() # Avoid text overlapping

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(pca.mydata, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")
#the function dimdesc() [in FactoMineR], for dimension description, can be used to identify the most significantly associated variables with a given principal component 
res.desc <- dimdesc(pca.mydata, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1
# Description of dimension 1
res.desc$Dim.2

# Biplot

pca.mydata3 <- PCA(df_to_cluster, quali.sup = 93:97, graph=FALSE)

fviz_pca_biplot(pca.mydata,
                col.ind = "#696969",  # Individuals color
                col.var = "contrib", 
repel = TRUE, alpha.var = "contrib") + scale_color_gradient2(low="white", mid="blue", 
                                                             high="red", midpoint=1) 

#col.ind = "#696969",  # Individuals color
#col.var = "#2E9FDF", # Variables color
#to color both individuals by clusters or k-means clusters and variables by contribution.
options(ggrepel.max.overlaps = 50)
p <- fviz_pca_biplot(pca.mydata,
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 5,
                label = "var",
                fill.ind = multiplex$DOS,
                col.ind = "white",
                # Color variable by groups
                col.var = "contrib", addEllipses = FALSE, invisible="quali", select.var = list(contrib = 50),
                repel = TRUE, alpha.var = "contrib", gradient.cols = magma(256, alpha = 1, begin = 0, end = 0.9, direction = 1))  +
  ggpubr::fill_palette(palette = c("#FDE725FF", "#f68f46ff", "#ad181c", "#593d9cff", "#3b2f1f")) + theme_light(base_size = 12)     # Indiviual fill color


p + geom_point(aes(shape = multiplex$Category), size=5) + scale_shape_manual(values=c(21, 22, 23)) #changes the shape accordingly to the groups 

ggpubr::show_point_shapes() 

#"#FDE725FF", "#55C667FF", "#238A8DFF","#453781FF", "black"


#ploting PCA only with the individuals data (not showing the  variables)
p2 <- fviz_pca_ind(pca.mydata,
                   # Fill individuals by groups
                   geom.ind = "point",
                   pointshape = 21,
                   pointsize = 6,
                   label = "ind",
                   fill.ind = multiplex$DOS, addEllipses =FALSE, invisible="quali",
                   col.ind = "white", col.var = "white",) + ggpubr::fill_palette(palette = c("#FDE725FF", "#55C667FF", "#238A8DFF","#453781FF", 'black')) + theme_light(base_size = 12) 

p2 + geom_point(aes(shape = multiplex$Category), size=6) + scale_shape_manual(values=c(21, 22, 23)) #changes the shape accordingly to the groups 
