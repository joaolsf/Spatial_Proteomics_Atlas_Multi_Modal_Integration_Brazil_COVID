library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(NbClust)
library(mclust)
library(M3C)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(viridis)

# reading file 
multiplex <- read.csv("COVID_Clinical_deceased_recovered2.csv")

# category 1
outcome <- multiplex$Outcome
cats <- ifelse(outcome <= 1, "Early death", 
                 ifelse(outcome > 1 & outcome <= 2, 'Late death', 'Recovered'))

# category 2
outcome2 <- multiplex$Outcome2
Disease_group <- ifelse(outcome2 <= 1, 'Fatal', 'Recovered')
               
# category 3
day <- multiplex$Day
Day_of_Symptoms <- ifelse(day<= 1, "Day 1", 
               ifelse(day > 1 & day <= 3, 'Day 2-6',
                      ifelse(day > 3 & day <= 5, 'Day 7',
                             ifelse(day > 5 & day <= 7, 'Day 8-13',
                                    ifelse(day > 7 & day <= 9, 'Day 14',
                                           ifelse(day > 9 & day <= 11, 'Day 15-27','Day 28'))))))

# retirando as colunas indesejadas para o heatmap
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("Age", "BMI", "Weight", "Category", "Category2", "Category3", "Outcome", "Outcome2", "Day", "DOS", "DOS2", "Disease_score")]

multiplex <- multiplex[,select_cols]
# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

# Prepare data
mydata <- scale(df_to_cluster)


#2- Elbow method
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:8) wss [i] <- sum(kmeans(mydata,centers=i)$withinss)
plot(1:8, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

#3- Model Based Clustering
library(mclust)
fit <- Mclust(mydata)
plot(fit) # plot results 
summary(fit) # display the best model

#4- Determine number of clusters using Monte Carlo reference-based consensus clustering (M3C)
#We recommend M3C only be used to cluster datasets with high numbers of samples (e.g. 60-1000). 
mydata_t <- t(mydata)
library(M3C)
pca(mydata_t,legendtextsize = 10,axistextsize = 10,dotsize=2)
res <- M3C(mydata_t, des = NULL, removeplots = TRUE, iters=20, clusteralg = 'km',
           objective='PAC', fsize=8, lthick=1, dotsize=1.25)

res$scores
res$plots[[1]]
res$plots[[2]]
res$plots[[4]]
res$plots[[3]]

data <- res$realdataresults[[3]]$ordered_data 
annon <- res$realdataresults[[3]]$ordered_annotation 
ccmatrix <- res$realdataresults[[3]]$consensus_matrix
head(annon)
pca(mydata_t,labels=annon$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2)

library(ComplexHeatmap)
ccl <- list()
x <- c("skyblue", "gold", "violet", "darkorchid", "slateblue", "forestgreen", 
                 "violetred", "orange", "midnightblue", "grey31", "black")
names(x) <- as.character(seq(1,11,by=1))
for (i in seq(2,10)){
    ccmatrix <- res$realdataresults[[i]]$consensus_matrix
    annon <- res$realdataresults[[i]]$ordered_annotation
n <- 10
seq <- rev(seq(0,255,by=255/(n)))
palRGB <- cbind(seq,seq,255)
mypal <- rgb(palRGB,maxColorValue=255)
ha = HeatmapAnnotation(
df= data.frame(Cluster=as.character(annon[,1])), col = list(Cluster=x))
ccl[[i]] <-  Heatmap(ccmatrix, name = "Consensus_index", top_annotation = ha, 
             col=mypal, show_row_dend = FALSE,
                show_column_dend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
                show_column_names = FALSE)
}
print(ccl[[i]])

pca(mydata_t,labels=annon$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2)

#5- K-Means Cluster Analysis
pdf("seeds_covid_fatal_recovered.pdf")
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
fit <- kmeans(mydata, 3) # 3 cluster solution
#plotting kmeans graph
library(cluster) 
library(fpc)
#plotcluster(mydata, fit$cluster)
clusplot(mydata, main = 3, fit$cluster, color=TRUE, shade=TRUE, labels=3, lines=0)
# get cluster means 
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata2 <- data.frame(mydata, fit$cluster)

# Fazendo o heatmap
library(ComplexHeatmap)
library(circlize)
## display a palettes simultanoeusly
#options of colors for heatmaps
library(RColorBrewer)
display.brewer.all(n=10, exact.n=FALSE)
col = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(3)
colorRamp2(c(-0.5, 0, 0.5), c("blue", "yellow", "red"), space = "RGB"),

cmeth <- c('ward.D', 'ward.D2', "complete", "single", "average", 
           "mcquitty", 'median', "centroid")
dmeth <- c("euclidean", "manhattan", "minkowski", "pearson")

#heatmap virus patterns/macrophage patterns/lymphocyte patterns

rannot <- rowAnnotation(df = data.frame(cats), 
                        col = list(cats = c("Early death" = "salmon", "Late death" = "light blue", "Recovered" = "grey")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))

rannot2 <- rowAnnotation(df = data.frame(Disease_group), 
                        col = list(Disease_group = c("Fatal" = "dark blue", "Recovered" = "grey")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))

rannot3 <- rowAnnotation(df = data.frame(Day_of_Symptoms), 
                        col = list(Day_of_Symptoms = c("Day 1" = "yellow","Day 2-6" = "orange", "Day 7" = "salmon", "Day 8-13" = "red", "Day 14" = "purple", "Day 15-27" = "black", "Day 28" = "brown")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))
#set.seed(1) 
#fit <- kmeans(mydata, 3) # 3 cluster solution

heat <- Heatmap(mydata, cluster_columns = TRUE, cluster_rows = FALSE,
                col = colorRamp2(c(1, 0.6, 0, -0.6, -1),brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                heatmap_legend_param = list(color_bar = "continuous"),
                show_row_dend = TRUE, show_column_dend = TRUE,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 14))
print(heat+rannot2+rannot3)

#1, 0.6, 0, -0.6, -1
#-1, -0.5, -0.25, 0, 0.25, 0.5, 1
#brewer.pal(n = 5, name = "Spectral"), space = "RGB")

 
pdf('heatmap4.pdf')
for (i in 1:40){
  heat <- Heatmap(mydata, clustering_distance_columns = "euclidean",
                  clustering_method_columns = "ward.D", km = 3, row_km_repeats = i,
                  col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), brewer.pal(n = 5, name = "YlOrRd"), space = "RGB"),
                  heatmap_legend_param = list(color_bar = "continuous"),
                  show_row_dend = TRUE, show_column_dend = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8))
  print(heat)
}
dev.off()


#Correlation deceased only
mydata2 <- mydata[1:323,]
#correlation only early death
mydata3 <- mydata[1:87,]
#correlation only late death
mydata4 <- mydata[88:323,]
#Correlation recovered only
mydata5 <- mydata[324:573,]


library(ggcorrplot)
corr <- round(cor(mydata2, use = "pairwise.complete.obs"), 2)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(corr)

#Visualize the correlation matrix
#method = "square" (default)
ggcorrplot(corr)

# Reordering the correlation matrix
# --------------------------------
# using hierarchical clustering
ggcorrplot(corr, hc.order = TRUE, outline.col = "grey",ggtheme = ggplot2::theme_gray, type = "lower", colors = c("#6D9EC1", "white", "#E46726"), tl.cex = 12,tl.srt = 90, lab = FALSE, p.mat = p.mat, sig.level = 0.05, insig = "blank")

# Change colors and theme
# --------------------------------
# Argument colors
ggcorrplot(corr, hc.order = TRUE,
           outline.col = "grey",
           ggtheme = ggplot2::theme_gray, tl.cex = 8,tl.srt = 90,
           colors = c("#6D9EC1", "white", "#E46726"), lab = TRUE, type = "lower")

library(Hmisc)
rcorr(mydata)

cor <- rcorr(mydata4, type = "spearman") # rcorr Calcula o p-value das correlacoes
cor$P
# plot das correlacoes de acordo com o p-value
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))

library(corrplot)
cor_matrix2 <- corrplot(cor$r, method = "square", p.mat = cor$P, sig.level = 0.05, 
                        insig = "blank", type = "full", is.corr=T, 
                        tl.cex=0.5, tl.col = "black", number.cex=0.5, number.font=0.5, diag=T, 
                        mar=c(1,1,1,1), order = "hclust", hclust.method = "ward.D", col = col4(200))



