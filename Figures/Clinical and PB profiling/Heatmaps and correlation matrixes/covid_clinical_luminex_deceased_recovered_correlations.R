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
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_correlation_longitudinal_samples_v2.csv")
#Reorder columns by biomarkers in the PB signatures
library(tibble)
multiplex <- as_tibble(multiplex)

col_order <- c("RecordID","PD_L1", "CXCL10", "CCL1", "MIP3b", "TF", "M_CSF", "TNFa", "MIP3a", "IL1RA", "IL27", "VCAM1", "MCP2", "MIP1b",
               "IL33", "IL18", "IL12p70", "IFNlambda2", "Ang2", "MCP1", "CXCL9", "Syndecan1", "OPN", "IL6", "ANG2_ANG1",
               "HGF", "vWFA2", "D_Dimer", "CXC3CL1", "MIP1a", "MPO", "IL8", 
               "IL11", "E_selectin", "IL1a", "TPO", "ICAM1",
               "IL1b", "IFNg", "IFNa", "IFNb", "IFNlambda3", "IL17", "IL21", "GranzymeB", "ADAMTS13", "FGFb", 
               "CXCL2", "CD40L", "EGF", "PDGF_BB", "Ang1", "PDGF_AA", "CXCL11", "RANTES", "CXCL1", "CXCL7", "C9",
               "SCF", "IL4", "IL7", "IL15", "IL5", "IL13", "G_CSF", "GM_CSF", "IL10", 
               "CCL11", "L_selectin", "CCL24", "CXCL4", "FasL", "Fibronectin",
               "Age",	"Weight", "BMI", "Hemoglobin", "Hematocrit", "Lymphocytes",	"Neutrophils", "NCLR", "Platelets",	"ALT", "Glucose",
               "Creatinine", "Urea", "LDH", "CRP", "Temperature", "HR", "RR", "SpO2", "DOSAA", "DOSAPA", "DDOUO", "DICU",	"DMV", "DOH",	
               "Disease_Score",	"Category",	"Category2", "Category3",	"Day")

multiplex <- multiplex[, col_order]
multiplex
write.csv(multiplex , "COVID_Clinical_Luminex_deceased_recovered_correlation_longitudinal_samples.csv")

# reading file (importar o arquivo)
multiplex <- read.csv("COVID_Clinical_Luminex_deceased_recovered_correlation_longitudinal_samples_v2.csv")

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
select_cols <- select_cols[!select_cols %in% c("Age", "BMI", "Weight",  "DOSAA","DOSAPA", "Disease_Score","Category", "Category2", "Category3", "Day")]

multiplex <- multiplex[,select_cols]
# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

# Prepare data
mydata <- scale(df_to_cluster)

#correlation only early death
#mydata_ed <- mydata[1:11,]
#correlation only late death
#mydata_ld <- mydata[12:47,]
#Correlation recovered only
mydata_recov <- mydata[64:123,]
#Correlation fatal 
mydata_fatal <- mydata[1:63,]

library(Hmisc)
cor <- rcorr(mydata_fatal, type = "pearson") # rcorr Calcula o p-value das correlacoes

library(tabletools)
#calculate adjusted p-values of the correlations: https://rdrr.io/github/JMLuther/tabletools/man/rcorr_padjust.html
cor2 <- rcorr_padjust(cor) # BH by default

#Extract r and p values to plot a corr plot for the correlation of specific variables
cor$r 
cor2$r
write.csv(cor$r , "Correlation_person_rcorr_clinical_luminex_fatal_longitudinal_samples_rvalues_v2.csv")

cor$P #unadjusted p-values
cor2$P #adjusted p-values using the BH correction method
write.csv(cor$P , "Correlation_pearson_rcorr_clinical_luminex_fatal_longitudinal_samples_pvalues_v2.csv")
write.csv(cor2$P , "Correlation_pearson_rcorr_clinical_luminex_fatal_longitudinal_samples_adjustpvalues_v2.csv")

# reading file (importar o arquivo)
multiplex2 <- read.csv("Correlation_person_rcorr_clinical_luminex_fatal_longitudinal_samples_rvalues_v2.csv")
multiplex3 <- read.csv("Correlation_pearson_rcorr_clinical_luminex_fatal_longitudinal_samples_adjustpvalues_v2.csv")

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

