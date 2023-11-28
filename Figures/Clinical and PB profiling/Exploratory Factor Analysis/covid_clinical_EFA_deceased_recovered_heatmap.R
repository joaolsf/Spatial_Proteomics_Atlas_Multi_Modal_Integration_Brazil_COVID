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
multiplex <- read.csv("COVID_Clinical_deceased_recovered_factor_loadings_iteration1.csv")

# data frame 
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$Clinical

# Prepare data
mydata <- (df_to_cluster)

# Heatmap
## display a palettes simultanoeusly
#options of colors for heatmaps
display.brewer.all(n=10, exact.n=FALSE)
col = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(3)
colorRamp2(c(-0.5, 0, 0.5), c("blue", "yellow", "red"), space = "RGB")

heat <- Heatmap(mydata, cluster_columns = FALSE, cluster_rows = FALSE,
                col = colorRamp2(c(0.7, 0.5, 0.3, 0, -0.3), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                heatmap_legend_param = list(color_bar = "continuous"),
                show_row_dend = TRUE, show_column_dend = TRUE,
                row_names_gp = gpar(fontsize = 14),
                column_names_gp = gpar(fontsize = 12), column_names_rot = 0)

print(heat)



