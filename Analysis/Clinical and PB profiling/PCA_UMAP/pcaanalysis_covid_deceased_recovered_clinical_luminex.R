# reading file (importar o arquivo)
multiplex <- read.csv("COVID_Luminex_fatal_recovered_PCA.csv")

# data frame para cluster
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

#quality and contribution individuals
#it’s also possible to color individuals by their cos2 values:
fviz_pca_ind(pca.mydata, col.ind = "cos2", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE, axes.linetype = "blank") + theme_bw(base_size = 14)
#You can also show the contrib of the corresponding individuals:
fviz_pca_ind(pca.mydata, col.ind = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, axes.linetype = "blank") + theme_bw(base_size = 14)
#You can also change the point size according the cos2 and colour based on contribution of the corresponding individuals:
fviz_pca_ind(pca.mydata, geom.ind = "point", col.ind = "contrib", pointsize = "cos2", 
             pointshape = 19, gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, axes.linetype = "blank") + theme_bw(base_size = 14)

# how to color individuals by groups
fviz_pca_ind(pca.mydata,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cats, # color by groups
             palette = c("black","grey", "#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")
#If you want confidence ellipses instead of concentration ellipses, use ellipse.type = “confidence”.
p <- fviz_pca_ind(pca.mydata,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = cats3, # color by groups
             palette = c("orange", "black","grey", "#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             legend.title = "Groups", pointshape = 19, mean.point = FALSE, axes.linetype = "blank")

p + geom_point(aes(colour = factor(cats3)), size = 5) + ggtitle("PC Patients") + theme_bw(base_size = 14) 

#paletes examples
c("black","#00AFBB", "#E7B800", "#FC4E07")
"npg"
"ggsci"
#Size and shape of plot elements
#Type ggpubr::show_point_shapes() to see available point shapes.
ggpubr::show_point_shapes() 
#The argument axes.linetype can be used to specify the line type of axes. Default is “dashed”. 
#Allowed values include “blank”, “solid”, “dotted”, etc. To see all possible values type ggpubr::show_line_types() in R.
ggpubr::show_line_types()

#To specify supplementary individuals and variables and qualitative variables, the function PCA() can be used as follow:
pca.mydata2 <- PCA(df_to_cluster, quali.sup = 33, graph=FALSE)

p <- fviz_pca_ind(pca.mydata2, col.ind.sup = "blue", repel = TRUE)
p <- fviz_add(p, pca.mydata2$quali.sup$coord, color = "red")
p
#To plot the graph with 3 clusters
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 33, pointsize = 5, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.65,
             palette = c("grey", "brown", "black"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with H&E patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 31, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("grey", "brown", "black", "blue","red"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with Spike lung patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 32, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("grey", "brown", "black"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with ISH patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 33, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("grey", "brown", "black", "blue"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with MAC387 patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 34, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("blue", "brown", "black", "grey"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with CD3 patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 35, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("blue", "brown", "black", "grey"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with CD20 patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 36, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("blue", "brown", "black", "grey", "red"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with Bronchial patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 37, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("grey", "brown", "black"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with Alveolar patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 38, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("grey", "brown", "black"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with Interstitial patterns
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 39, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("grey", "brown", "black"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)

#To plot the graph with Macrophage clusters
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 40, pointsize = 10, pointshape = 19,
             addEllipses =FALSE, ellipse.type = "confidence", ellipse.level = 0.95,
             palette = c("red", "brown", "black", "blue", "grey"), repel = TRUE, mean.point = FALSE, axes.linetype = "blank") + theme_bw(base_size = 14)





"#00BA38", "#619CFF", "#F8766D"

#All the outputs of the PCA (individuals/variables coordinates, contributions, etc) can be exported at once, into a TXT/CSV file, using the function write.infile() [in FactoMineR] package:
# Export into a TXT file
write.infile(pca.mydata, "pca.txt", sep = "\t")
# Export into a CSV file
write.infile(pca.mydata, "pca.csv", sep = ";")



#Biplot

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



#c("salmon", "blue", "dark grey")

#+ scale_color_gradient2(low="orange", mid="purple", 
                        high="black", midpoint=5) +
#colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(7))
#gradient.cols = magma(5, direction = -1)
#viridis
#magma
#plasma
#inferno
#mako

#brewer.pal(n = 5, name = "YlOrRd"), space = "RGB"
#palette = c("salmon", "light blue")
#"grey"
#"dark green", 
#"black"
#"blue", "orange","dark green","red"
#"grey", "red", "brown","light blue","blue", "dark blue"
