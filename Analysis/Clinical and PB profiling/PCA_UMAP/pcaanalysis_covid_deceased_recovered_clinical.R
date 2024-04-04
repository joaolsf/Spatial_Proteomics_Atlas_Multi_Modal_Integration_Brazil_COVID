
multiplex <- read.csv("COVID_Clinical_deceased_recovered.csv")
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

#subsetting active variables for the PCA
df_to_cluster.active <- df_to_cluster[,1:22]

# Prepare data
mydata <- scale(df_to_cluster.active)

#1- PCA Analysis
library("FactoMineR")
library("factoextra")
library(RColorBrewer)
library(ggrepel)
library(viridis)
library(ggplot2)
library(umap)
library(ggthemes)

#1- PCA with prcomp function
pca <- as.data.frame(prcomp(mydata)$x) 
plot(pca$PC1, pca$PC2)

#1.1- K-Means Cluster Analysis
pdf("seeds_kmeans_covid_clinica_prcomp.pdf")
for(i in 1:40){
  set.seed(i)
  fit <- kmeans(mydata, 3) # 3 cluster solution
  #plotting kmeans graph
  library(cluster) 
  library(fpc)
  #plotcluster(mydata, fit$cluster)
  clusplot(mydata, main = i, fit$cluster, color=TRUE, shade=TRUE, labels=3, lines=0)
}
dev.off()

set.seed(1) 
fit <- kmeans(mydata, 2) # 3 cluster solution
#plotting kmeans graph
library(cluster) 
library(fpc)
#plotcluster(mydata, fit$cluster)
clusplot(mydata, main = 2, fit$cluster, color=TRUE, shade=TRUE, labels=3, lines=0)
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata2 <- data.frame(mydata, fit$cluster)

#merging mydata e mydata2 para fazer o plot de PCA com os clusters
pca2 <- data.frame(subject = rownames(pca), PC1 = pca$PC1, PC2 = pca$PC2)
mydata3 <- data.frame(subject = rownames(mydata2), kmeans_clustering = mydata2$fit.cluster)
merged <- merge(pca2, mydata3, by = "subject", all = TRUE)
merged2 <- data.frame(multiplex, subject = rownames(multiplex), kmeans_clustering = merged$kmeans_clustering, PC1 = merged$PC1, PC2 = merged$PC2)

library(ggplot2)
library(ggthemes)
library(ggrepel)
ggplot(pca2, aes(x = PC1, y = PC2)) + geom_point()
p <- ggplot(pca2, aes(x = PC1, y = PC2, label = subject))
p + geom_point(aes(colour = factor(multiplex$Category2)), size = 5) + ggtitle("PC Patients") + theme_bw(base_size = 20) 

ggplot(merged, aes(x = PC1, y = PC2)) + geom_point()
p <- ggplot(merged, aes(PC1, PC2, label = subject, col=factor(kmeans_clustering)))

myColors <- c("grey", "salmon","dark blue") 
p + geom_point (aes(shape = factor(multiplex$Category2)), size = 3) + scale_color_manual(values=myColors) + ggtitle("PC Patients") + theme_bw(base_size = 20) 

#2- PCA with PCA function
pca.mydata <- PCA(mydata, scale.unit = FALSE, ncp = 10, graph = TRUE)
print(pca.mydata)

#2.1- extract the eigenvalues/variances of PCs
eig.val <- get_eigenvalue(pca.mydata)
eig.val
#visualize the eigenvalues
fviz_eig(pca.mydata)

#2.2- Visualize the results individuals and variables, respectively.
fviz_pca_ind(pca.mydata)
fviz_pca_var(pca.mydata)
fviz_eig(pca.mydata, addlabels = TRUE, ylim = c(0, 50))

#2.3- Method to extract the results, for variables,
get_pca_var(pca.mydata)
var <- get_pca_var(pca.mydata)
var
fviz_pca_var(pca.mydata, col.var = "black")

#2.4- Method to extract the results for individuals
get_pca_ind(pca.mydata)
ind <- get_pca_ind(pca.mydata)
ind
fviz_pca_ind(pca.mydata, col.var = "black")

#2.5- The quality of representation of the variables on factor map is called cos2 (square cosine, squared coordinates)
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

#2.6- The function dimdesc() [in FactoMineR], for dimension description, can be used to identify the most significantly associated variables with a given principal component 
res.desc <- dimdesc(pca.mydata, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1
# Description of dimension 2
res.desc$Dim.2

write.csv(res.desc$Dim.1, "PC1_clinical_correlation.csv")
write.csv(res.desc$Dim.2, "PC2_clinical_correlation.csv")

#plot factor loadings for PC1
weights <- read.csv("PC1_clinical_correlation.csv")
plot_ <- ggplot(weights,
                aes(x= reorder(Feature,
                               (Correlation), decreasing = TRUE), y = (Correlation), fill = Category)) +
  geom_bar(stat = "identity", width=0.7, position=position_dodge(width=1)) +
  theme_classic(base_size = 18) +
  xlab("") + 
  ylab("PC1 Loadings") + scale_fill_manual("legend", values = c("Disease severity" = "black", "Tissue injury marker" = "#F79C79", "Hematological parameter" = "#800080", "Phisiology" = 'salmon', "Acute inflammatory marker" = '#21908CFF')) +
  theme(plot.title = element_text(size=18), 
        axis.text.y = element_text(face = "plain", size=18, angle=90, vjust=.5, hjust=1),
        axis.text.x = element_text(face = "plain", size=22, angle=90, vjust=.5, hjust=1))
plot_

#plot factor loadings for PC2
weights2 <- read.csv("PC2_clinical_correlation.csv")
plot_ <- ggplot(weights2,
                aes(x= reorder(Feature,
                               (Correlation), decreasing = TRUE), y = (Correlation), fill = Category)) +
  geom_bar(stat = "identity", width=0.7, position=position_dodge(width=1)) +
  theme_classic(base_size = 18) +
  xlab("") + 
  ylab("PC2 Loadings") + scale_fill_manual("legend", values = c("Disease severity" = "black", "Tissue injury marker" = "#F79C79", "Hematological parameter" = "#800080", "Phisiology" = 'salmon', "Acute inflammatory marker" = '#21908CFF')) +
  theme(plot.title = element_text(size=18), 
        axis.text.y = element_text(face = "plain", size=18, angle=90, vjust=.5, hjust=1),
        axis.text.x = element_text(face = "plain", size=22, angle=90, vjust=.5, hjust=1))
plot_

#2.7- KMEANS CLUSTERING in the PCA REDUCTION
# Create a grouping of individuals using kmeans
# Create 3 groups of individuals (centers = 3)
set.seed(123)
res.km <- kmeans(ind$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_ind(pca.mydata, col.ind = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")

res.km$cluster
multiplex2 <- data.frame(multiplex, subject = rownames(multiplex), kmeans_cluster = res.km$cluster)
write.csv(multiplex2, 'COVID_Clinical_deceased_recovered_kmeans.csv')

#2.8- Quality and contribution individuals
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

#2.9- How to color individuals by groups
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

#Size and shape of plot elements
#Type ggpubr::show_point_shapes() to see available point shapes.
ggpubr::show_point_shapes() 
#The argument axes.linetype can be used to specify the line type of axes. Default is “dashed”. 
#Allowed values include “blank”, “solid”, “dotted”, etc. To see all possible values type ggpubr::show_line_types() in R.
ggpubr::show_line_types()

#To specify supplementary individuals and variables and qualitative variables, the function PCA() can be used as follow:
pca.mydata2 <- PCA(df_to_cluster, quali.sup = 21:28, graph=FALSE)

p <- fviz_pca_ind(pca.mydata2, col.ind.sup = "blue", repel = TRUE)
p <- fviz_add(p, pca.mydata2$quali.sup$coord, color = "red")
p
#To plot the graph with 3 clusters
fviz_pca_ind(pca.mydata2, geom.ind = "point", habillage = 21, pointsize = 5, pointshape = 19,
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

# Export into a TXT file
write.infile(pca.mydata, "pca.txt", sep = "\t")
# Export into a CSV file
write.infile(pca.mydata, "pca.csv", sep = ";")

#2.10- PCA Biplot
pca.mydata3 <- PCA(df_to_cluster, quali.sup = 27:31, graph=FALSE)

fviz_pca_biplot(pca.mydata,
                col.ind = "#696969",  # Individuals color
                col.var = "contrib", 
repel = TRUE, alpha.var = "contrib") + scale_color_gradient2(low="white", mid="blue", 
                                                             high="red", midpoint=1) 
#color individuals by disease outcome.
p <- fviz_pca_biplot(pca.mydata,
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 4,
                label = "var",
                fill.ind = multiplex$Category2,
                col.ind = "white",
                # Color variable by groups
                col.var = "contrib", addEllipses = FALSE, invisible="quali",
                repel = TRUE, alpha.var = "contrib", gradient.cols = magma(256, direction = 1, begin = 0, end = 0.85)) +
  ggpubr::fill_palette(palette = c("dark blue", "light grey")) + theme_bw(base_size = 24) 
p

# color individuals by disease progression.
p2 <- fviz_pca_biplot(pca.mydata,
                     # Fill individuals by groups
                     geom.ind = "point",
                     pointshape = 21,
                     pointsize = 4,
                     label = "var",
                     fill.ind = multiplex$Category,
                     col.ind = "white",
                     # Color variable by groups
                     col.var = "contrib", addEllipses = FALSE, invisible="quali",
                     repel = TRUE, alpha.var = "contrib", gradient.cols = magma(256, direction = 1, begin = 0, end = 0.85)) +
  ggpubr::fill_palette(palette = c("salmon", "light blue", "light grey")) + theme_bw(base_size = 24) 
p2

multiplex$kmeans_cluster <- as.factor(multiplex$kmeans_cluster)
#to colot individuals by k-means cluster and shape by disease group
p3 <- fviz_pca_biplot(pca.mydata,
                     # Fill individuals by groups
                     geom.ind = "point",
                     pointshape = 21,
                     pointsize = 4,
                     label = "var",
                     fill.ind = multiplex$kmeans_cluster,
                     col.ind = "white",
                     # Color variable by groups
                     col.var = "contrib", addEllipses = FALSE, invisible="quali",
                     repel = TRUE, alpha.var = "contrib", gradient.cols = magma(256, direction = 1, begin = 0, end = 0.85)) +
  ggpubr::fill_palette(palette = c("#EFC000FF","#0073C2FF", "#868686FF")) + theme_bw(base_size = 12) #color palette from kmeans function above
p3

p3 + geom_point(aes(shape = factor(multiplex$Category2), size=4)) + scale_shape_manual(values=c(21, 23)) +
  theme(axis.text.y = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5))

#another way to plot the graph above changing the shape accordingly to the groups - but it does not plot the variables on top
pca <- as.data.frame(prcomp(mydata)$x) 
pca2 <- data.frame(subject = rownames(pca), PC1 = pca$PC1, PC2 = pca$PC2)
library(ggplot2)
library(ggthemes)
library(ggrepel)
ggplot(pca2, aes(x = PC1, y = PC2)) + 
       geom_point(aes(shape=multiplex$Category2, fill=multiplex$DOS2, size=0.5)) + 
       scale_shape_manual(values=c(21, 22)) + 
       ggpubr::fill_palette(palette = c("#FDE725FF", "#DCE319FF", "#55C667FF", "#238A8DFF","#453781FF","#440154FF", "black")) + 
       theme_bw(base_size = 12)

#ploting PCA only with the individuals data (not showing the  variables)
p2 <- fviz_pca_ind(pca.mydata,
                     # Fill individuals by groups
                     geom.ind = "point",
                     pointshape = 21,
                     pointsize = 4,
                     label = "ind",
                     fill.ind = multiplex$DOS2, addEllipses =FALSE, invisible="quali",
                     col.ind = "white", col.var = "white",) + ggpubr::fill_palette(palette = c("#FDE725FF", "#DCE319FF", "#55C667FF", "#238A8DFF","#453781FF","#440154FF", "black")) + theme_light(base_size = 12) 

p2 + geom_point(aes(shape = multiplex$Category2), size=4) + scale_shape_manual(values=c(21, 23)) #changes the shape accordingly to the groups 

#3- UMAP
set.seed(1234)
mydata.umap <- umap(mydata) 
mydata.umap
head(mydata.umap$layout, 3)
mydata.labels <- multiplex[, "Category2"]


library(tidyverse)
umap_df <- mydata.umap$layout %>%
  as.data.frame() %>%
  rename(UMAP1="V1",
         UMAP2="V2")

umap_df %>% head()

#Plot UMAP by disease group
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = multiplex$Category2,
             shape = multiplex$Category2,
             size=3)) +
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") + scale_color_manual(values = c("dark blue", "grey")) + theme_light(base_size = 24)


  theme(legend.justification=c(),
        legend.position='bottom', plot.title = element_text(size=20), 
        axis.text.y = element_text(face = "plain", size=22, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=22, angle=0, vjust=1, hjust=0.5), 
        strip.text.x = element_text(size = 22, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank"))

#Plot UMAP by kmeans clusters
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = grp,
             shape = multiplex$Category2,
             size=3)) +
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") + scale_color_manual(values = c("#EFC000FF","#0073C2FF", "#868686FF")) + theme_light(base_size = 12)

#Plot UMAP by early vs late death
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = multiplex$Category,
             shape = multiplex$Category2,
             size=3)) +
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot") + scale_color_manual(values = c("salmon", "light blue", "grey")) + theme_light(base_size = 12) 


