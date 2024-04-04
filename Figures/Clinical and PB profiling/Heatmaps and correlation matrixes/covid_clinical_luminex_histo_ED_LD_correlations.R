library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(viridis)
library(longitudinal)
library(ggcorrplot)
library(Hmisc)
library(corrplot)

# reading file
multiplex <- read.csv("COVID_Clinical_PB_histology_factors_longitudinal_corr.csv")

# retirando as colunas indesejadas para o heatmap
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("BlockID")]
multiplex <- multiplex[,select_cols]

# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

# Prepare data
mydata <- scale(df_to_cluster)

library(Hmisc)
cor <- rcorr(mydata, type = "spearman") # rcorr Calcula o p-value das correlacoes

library(tabletools)
#calculate adjusted p-values of the correlations: https://rdrr.io/github/JMLuther/tabletools/man/rcorr_padjust.html
cor2 <- rcorr_padjust(cor) # BH by default

#Extract r and p values to plot a corr plot for the correlation of specific variables
cor$r 
cor2$r
write.csv(cor$r , "Correlation_spearman_rcorr_clinical_PB_histo_signatures_longitudinal_rvalues.csv")

cor$P #unadjusted p-values
cor2$P #adjusted p-values using the BH correction method
write.csv(cor$P , "Correlation_spearman_rcorr_clinical_PB_histo_signatures_longitudinal_pvalues.csv")
write.csv(cor2$P , "Correlation_spearman_rcorr_clinical_PB_histo_signatures_longitudinal_adjustpvalues.csv")

# reading file (importar o arquivo)
multiplex2 <- read.csv("Correlation_spearman_rcorr_clinical_PB_histo_signatures_longitudinal_rvalues.csv")
multiplex3 <- read.csv("Correlation_spearman_rcorr_clinical_PB_histo_signatures_longitudinal_pvalues.csv")

rownames(multiplex2) <- multiplex2$Clinical
rownames(multiplex3) <- multiplex3$Clinical

multiplex2 <- multiplex2[,-1]
multiplex2 <- as.matrix.data.frame(multiplex2)

multiplex3 <- multiplex3[,-1]
multiplex3 <- as.matrix.data.frame(multiplex3)

#multiplex4 <- t(multiplex2)
#multiplex5 <- t(multiplex3)

#plot the correlation of specific variables using the files above
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))
library(corrplot)
cor_matrix <- corrplot(multiplex2, method = "square", p.mat = multiplex3, 
                       type = "full", is.corr=T, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.75,pch.col = 'white',
                       tl.cex=0.8, tl.col = "black", number.cex=0.5, number.font=0.5, diag=F,cl.cex = 1, cl.ratio = 0.1,
                       mar=c(0.05,0.05,0.05,0.05), col = col4(200), addgrid.col="light grey")


#order = "hclust", hclust.method = "ward.D"
#insig = 'blank',


library("longitudinal")
mydata2 <- as.longitudinal(mydata2, repeats=c(11,7,10,9,7,3), time=c(1,3,5,7,11,14))
mydata5 <- as.longitudinal(mydata5, repeats=c(10,10,10,10,10,10), time=c(1,3,5,7,11,14))

dynpc <- dyn.cor(mydata2, lambda = 0.4798)
dynpc2 <- as.data.frame(dynpc)

