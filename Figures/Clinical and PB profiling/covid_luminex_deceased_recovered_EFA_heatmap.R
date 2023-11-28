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
multiplex <- read.csv("COVID_Luminex_deceased_recovered_factor_loadings.csv")

# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$CytokineID

# Prepare data
mydata <- df_to_cluster
mydata2 <- as.matrix(mydata)

factor1 <- mydata[,1:10]
factor2 <- mydata[,11:19]
factor3 <- mydata[,20:30]
factor4 <- mydata[,31:61]
factor5 <- mydata[,62:66]
  

# Fazendo o heatmap
library(ComplexHeatmap)
library(circlize)
## display a palettes simultanoeusly
#options of colors for heatmaps
library(RColorBrewer)
display.brewer.all(n=10, exact.n=FALSE)
col = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(3)
colorRamp2(c(-0.5, 0, 0.5), c("blue", "yellow", "red"), space = "RGB")

rannot <- rowAnnotation(df = data.frame(cats), 
                        col = list(cats = c("Early death" = "salmon", "Late death" = "light blue", "Recovered" = "grey")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))

rannot2 <- rowAnnotation(df = data.frame(cats2), 
                        col = list(cats2 = c("Day 1" = "yellow","Day 3" = "orange", "Day 5" = "salmon", "Day 7" = "red", "Day 11" = "purple", "Day 14" = "black")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))

heat <- Heatmap(mydata2, cluster_columns = FALSE, cluster_rows = FALSE,
                col = colorRamp2(c(0.7, 0.5, 0.3, 0, -0.3), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                heatmap_legend_param = list(color_bar = "continuous"),
                show_row_dend = TRUE, show_column_dend = TRUE,
                row_names_gp = gpar(fontsize = 12),
                column_names_gp = gpar(fontsize = 12), column_names_rot = 360, column_names_centered = TRUE)

print(heat)

