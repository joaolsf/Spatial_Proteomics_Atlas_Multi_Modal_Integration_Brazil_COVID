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

# reading file (importar o arquivo)
multiplex <- read.csv("COVID_Clinical_PB_histology_factors_only_last_sample_corr.csv")
#multiplex <- multiplex[,1:8]
# category 1
outcome <- multiplex$Category2
cats <- ifelse(outcome <= 1, "Early death", 
                 ifelse(outcome > 1 & outcome <= 2, 'Late death', 'Recovered'))
               
# category 2
day <- multiplex$Category3
cats2 <- ifelse(day<= 1, "Day 1", 
               ifelse(day > 1 & day <= 3, 'Day 3',
                      ifelse(day > 3 & day <= 5, 'Day 5',
                             ifelse(day > 5 & day <= 7, 'Day 7',
                                    ifelse(day > 7 & day <= 11, 'Day 11', 'Day 14')))))

# retirando as colunas indesejadas para o heatmap
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("Age", "BMI", "Weight", "Category", "Category2", "Category3", "Outcome", "Day",
                                               "Hemoglobin", "Hematocrit", "Lymphocytes", "Neutrophils","NCLR", "Platelets",
                                               "ALT", "Glucose", "Creatinine", "Urea", "LDH", "CRP", "Temperature", "HR", "RR", "SpO2", "DOSAA", "DOSAPA", "DDOUD", "DICU", "DMV", "DOH", "Disease_Score", "DOS")]

multiplex <- multiplex[,select_cols]
# data frame para cluster
df_to_cluster <- multiplex[,-c(1:2)]
rownames(df_to_cluster) <- multiplex$RecordID

# Prepare data
mydata <- scale(df_to_cluster)
#mydata <- mydata[,1:5]

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

#subsetting mydata accordingly to day of sampling within each disease group
#late death
#day1
mydata2 <- mydata[12:18,]
#day3
mydata3 <- mydata[19:24,]
#day5
mydata4 <- mydata[25:31,]
#day7
mydata5 <- mydata[32:38,]
#day11
mydata6 <- mydata[39:44,]

#recovered
#day1
mydata7 <- mydata[48:57,]
#day3
mydata8 <- mydata[58:67,]
#day5
mydata9 <- mydata[68:77,]
#day7
mydata10 <- mydata[78:87,]
#day11
mydata11 <- mydata[88:97,]
#day14
mydata12 <- mydata[98:107,]

rannot <- rowAnnotation(df = data.frame(cats), 
                        col = list(cats = c("Early death" = "salmon", "Late death" = "light blue", "Recovered" = "grey")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))

rannot2 <- rowAnnotation(df = data.frame(cats2), 
                        col = list(cats2 = c("Day 1" = "yellow","Day 3" = "orange", "Day 5" = "salmon", "Day 7" = "red", "Day 11" = "purple", "Day 14" = "black")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))

heat <- Heatmap(mydata, cluster_columns = TRUE, cluster_rows = FALSE,clustering_distance_columns = "euclidean",
                clustering_method_columns = "ward.D",
                col = colorRamp2(c(1, 0.6, 0, -0.6, -1), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                heatmap_legend_param = list(color_bar = "continuous"),
                show_row_dend = TRUE, show_column_dend = TRUE,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 10))

print(heat+rannot+rannot2)

#1, 0.5, 0.2, 0, -0.2,-0.5, -1
#brewer.pal(n = 5, name = "RdYlBu"), space = "RGB")
#brewer.pal(n = 5, name = "YlOrRd"), space = "RGB")
#magma(5, alpha = 1, begin = 0, end = 1, direction = -1)


set.seed(1) 
fit <- kmeans(mydata, 3) # 3 cluster solution
pdf('heatmap_covid_fatal_clusters3.pdf')
for(cme in cmeth){
  for(dme in dmeth){
    heat <- Heatmap(mydata, cluster_columns = TRUE, cluster_rows = FALSE, clustering_distance_columns = dmeth,
                    clustering_method_columns = cmeth,
                    name = paste0(cme, ' ', dme), 
                    col = colorRamp2(c(1, 0.6, 0, -0.6, -1), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                    heatmap_legend_param = list(color_bar = "continuous"),
                    show_row_dend = TRUE, show_column_dend = FALSE,
                    row_names_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 6))
    
    print(heat++rannot+rannot2)
  }
}
dev.off()

pdf('heatmap_covid.pdf')
for (i in 1:40){
    heat <- Heatmap(mydata, cluster_columns = TRUE, cluster_rows = FALSE, clustering_distance_columns = "euclidean",
                    clustering_method_columns = "ward.D",
                    name = paste0(cme, ' ', dme), 
                    col = colorRamp2(c(1, 0.6, 0, -0.6, -1), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                    heatmap_legend_param = list(color_bar = "continuous"),
                    show_row_dend = TRUE, show_column_dend = FALSE,
                    row_names_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 6))
    
    print(heat+rannot+rannot2)
  }
dev.off()
 
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

cor <- rcorr(mydata, type = "spearman") # rcorr Calcula o p-value das correlacoes
cor$P
# plot das correlacoes de acordo com o p-value
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))

library(corrplot)
cor_matrix2 <- corrplot(cor$r,method = "square", p.mat = cor$P, 
                        type = "full", is.corr=T, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,pch.col = 'white',
                        tl.cex=0.5, tl.col = "black", number.cex=0.5, number.font=0.5, diag=F,cl.cex = 0.5, cl.ratio = 0.1,
                        mar=c(0.05,0.05,0.05,0.05), col = col4(200), addgrid.col="light grey")

#order = "hclust", hclust.method = "ward.D",
#Extract r and p values to plot a corr plot for the correlation of specific variables
cor$r 
write.csv(cor$r , "Correlation_spearman_rcorr_clinical_PB_histo_signatures_only_last_sample_rvalues.csv")
cor$P
write.csv(cor$P , "Correlation_spearman_rcorr_clinical_PB_histo_signatures_only_last_sample_pvalues.csv")

# reading file (importar o arquivo)
multiplex2 <- read.csv("Correlation_spearman_rcorr_clinical_PB_histo_signatures_only_last_sample_rvalues.csv")
multiplex3 <- read.csv("Correlation_spearman_rcorr_clinical_PB_histo_signatures_only_last_sample_pvalues.csv")

rownames(multiplex2) <- multiplex2$X
rownames(multiplex3) <- multiplex3$X

multiplex2 <- multiplex2[,-1]
multiplex2 <- as.matrix.data.frame(multiplex2)

multiplex3 <- multiplex3[,-1]
multiplex3 <- as.matrix.data.frame(multiplex3)

cor_matrix3 <- corrplot(multiplex2,method = "square", p.mat = multiplex3, 
                        type = "full", is.corr=T, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 2,pch.col = 'white',
                        tl.cex=1, tl.col = "black", number.cex=2, number.font=2, diag=F,cl.cex = 0.5, cl.ratio = 0.1,
                        mar=c(0.05,0.05,0.05,0.05), col = col4(200), addgrid.col="light grey")



