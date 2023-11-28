# reading file (importar o arquivo)
multiplex <- read.csv("COVID_Clinical_Luminex_Histology_IMC_early_late_death.csv")

length <- multiplex$Length
cats <- ifelse(length <= 1, 'Early death', 'Late death')
# retirando as colunas indesejadas for correlation
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("BlockID", "Length","Score")]

#"MAC387_Bronchial_Pattern", "MAC387_Alveolar_Pattern", "MAC387_Interstitial_Pattern", "Angiogenesis_Intussusception", "Angiogenesis_Sprouting", "Vascular_Congestion", "Arterial_Thrombi", "Venous_Thrombi", "Interstitial_Oedema","Lymphocytes_HE", "Megakaryocytes_HE
multiplex <- multiplex[,select_cols]

# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

#subsetting active variables for the PCA
df_to_cluster.active <- df_to_cluster[,1:237]

# Prepare data
mydata <- (df_to_cluster.active)

#1-Handling missing data in the dataframe:
library("missMDA")
library("FactoMineR")
library("factoextra")
library(RColorBrewer)
library(ggrepel)
library(viridis)

#select the number of dimensions that will be used in the algorithm using the function estim_ncpPCA
#By default, the function estim_ncpPCA uses the GCV method.
#It is possible to use another cross-validation strategy by specifying the argument as follows:
#ncomp$ncp <- estim_ncpPCA(geno, ncp.min = 0, ncp.max = 6, method.cv="Kfold", nbsim=100,pNA=0.05)
#this piece of the code takes forever and it has not been working...it gives an error with lots of NA values.
#ncomp <- estim_ncpPCA(mydata)
#ncomp$ncp

#Perform the (regularized) iterative PCA algorithm with the number of dimensions selected in the previous step, using the function imputePCA:
res.imp <- imputePCA(mydata, ncp = 6)

#perform PCA on the imputed data set. To this end, we propose to use the PCA function of the FactoMineR package:
pca.mydata <- PCA(res.imp$completeObs)
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
fviz_contrib(pca.mydata, choice = "var", axes = 1, top = 50, xtickslab.rt = 90, fill = "#800080", color = "#800080", ggtheme = theme_light())
# Contributions of variables to PC2
fviz_contrib(pca.mydata, choice = "var", axes = 2, top = 50,  xtickslab.rt = 90, fill = "#800080", color = "#800080", ggtheme = theme_light())
#The total contribution to PC1 and PC2 is obtained with the following R code:
fviz_contrib(pca.mydata, choice = "var", axes = 1:5, top = 50,  xtickslab.rt = 90, fill = "#800080")
# Color by contrib values: quality on the factor map
fviz_pca_var(pca.mydata, col.var = "contrib", select.var = list(contrib = 30),
             repel = TRUE, alpha.var = "contrib") + scale_color_gradient2(low="white", mid="blue", 
            high="red", midpoint=0.9) + theme_light() # Avoid text overlapping

#All the outputs of the PCA (individuals/variables coordinates, contributions, etc) can be exported at once, into a TXT/CSV file, using the function write.infile() [in FactoMineR] package:
# Export into a TXT file
write.infile(pca.mydata, "pca.txt", sep = "\t")
# Export into a CSV file
write.infile(pca.mydata, "pca.csv", sep = ";")

#PCA Biplot
#to color both individuals by clusters or k-means clusters and variables by contribution.
options(ggrepel.max.overlaps = 50)
display.brewer.all(n=10, exact.n=FALSE)

fviz_pca_biplot(pca.mydata,
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 6,
                label = "var",
                fill.ind = multiplex$DOS_PCA,
                col.ind = "white",
                # Color variable by groups
                col.var = "contrib", select.var = list(contrib = 50), addEllipses =FALSE, invisible="quali",
                repel = TRUE, alpha.var = "contrib", gradient.cols = mako(256, alpha = 1, begin = 0, end = 0.95, direction = -1)) +
  ggpubr::fill_palette(palette = magma(5, alpha = 1, begin = 0, end = 0.95, direction = -1)) + theme_bw(base_size = 12)  # Indiviual fill color

magma(6, alpha = 1, begin = 0, end = 0.95, direction = -1)
+ scale_color_gradient2(low="orange", mid="purple", 
                                   high="black", midpoint=1.5) +
  
colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(7))

ggpubr::fill_palette(palette = c("salmon", "light blue"))

palette = magma(6, alpha = 1, begin = 0, end = 0.95, direction = -1)

"grey"
"dark green", 
"black"
"blue", "orange","dark green","red"
"grey", "red", "brown","light blue","blue", "dark blue"

p2 <- fviz_pca_ind(pca.mydata, 
             geom.ind = "point",
             pointshape = 21,
             pointsize = 5,
             label = "ind",
             fill.ind = multiplex$Day_PCA, addEllipses =FALSE, invisible="quali", ellipse.type = "confidence", ellipse.level = 0.95,
             col.ind = "white") + ggpubr::fill_palette(palette = magma(6, alpha = 1, begin = 0, end = 1, direction = -1)) + theme_bw(base_size = 12)

p2 + geom_point(aes(shape = multiplex$Category), size=5) + scale_shape_manual(values=c(21, 22)) #changes the shape accordingly to the groups 



