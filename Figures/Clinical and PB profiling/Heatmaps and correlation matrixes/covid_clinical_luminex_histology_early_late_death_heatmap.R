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
multiplex <- read.csv("Correlation_spearman_rcorr_clinical_PB_histo_signatures_only_last_sample_rvalues.csv")
multiplex2 <- read.csv("Correlation_spearman_rcorr_clinical_PB_histo_signatures_only_last_sample_pvalues.csv")

rownames(multiplex) <- multiplex$X
rownames(multiplex2) <- multiplex2$X

multiplex <- multiplex[,-1]
multiplex <- as.matrix.data.frame(multiplex)

multiplex2 <- multiplex2[,-1]
multiplex2 <- as.matrix.data.frame(multiplex2)

# Histopattern_lung
parasit <- multiplex$Histopattern_lung
cats <- ifelse(parasit >= 0 & parasit <= 1, 'ND',
               ifelse(parasit > 1 & parasit <= 2, 'Acute_severe_fibrin_high_neutros', 
                      ifelse(parasit > 2 & parasit <= 3, 'Subacute_mononuclear_inflammatory_cells', 'Haemorrhagic')))

# DAD phase
dad <- multiplex$DAD_phase
cats2 <- ifelse(dad >= 0 & dad < 1, 'ND',
                ifelse(dad >= 1 & dad < 2, 'Exudative',
                 ifelse(dad >= 2 & dad < 3, 'Proliferative', 'Fibrotic')))

# spike IHC lung
spike <- multiplex$IHC_Virus_Lung
cats3 <- ifelse(spike >= 0 & spike < 1, 'ND',
              ifelse(spike >= 1 & spike < 2, 'Negative', 'Positive'))

# virus ISH lung
ish <- multiplex$ISH_Virus_Lung
cats4 <- ifelse(ish >= 0 & ish < 1, 'ND',
                ifelse(ish >= 1 & ish < 2, 'Negative', 'Positive'))
                     
# Cellularity
mm <- multiplex$Cellularity
cats30 <- ifelse(mm >= 99 & mm < 100,'ND',
                 ifelse(mm >= 0 & mm < 1,'Normal',
                        ifelse(mm >=1 & mm < 2,'Mild',
                               ifelse(mm >= 2 & mm < 3, 'Moderate','High'))))

# Monocytes_Macrophages_HE
mm <- multiplex$Monocytes_Macrophages_HE
cats31 <- ifelse(mm >= 99 & mm < 100,'ND',
                 ifelse(mm >= 0 & mm < 1,'Normal',
                        ifelse(mm >=1 & mm < 2,'Mild',
                               ifelse(mm >= 2 & mm < 3, 'Moderate','High'))))

# Neutrophils_HE
neutro <- multiplex$Neutrophils_HE
cats32 <- ifelse(neutro >= 99 & neutro < 100,'ND',
                 ifelse(neutro >= 0 & neutro < 1,'Normal',
                        ifelse(neutro >=1 & neutro < 2,'Mild',
                               ifelse(neutro >= 2 & neutro < 3, 'Moderate','High'))))

# Lymphocytes_HE
lt <- multiplex$Lymphocytes_HE
cats33 <- ifelse(lt >= 99 & lt < 100,'ND',
                 ifelse(lt >= 0 & lt < 1,'Normal',
                        ifelse(lt >=1 & lt < 2,'Mild',
                               ifelse(lt >= 2 & lt < 3, 'Moderate','High'))))

# Plasma_cells
pc <- multiplex$Plasma_cells
cats34 <- ifelse(pc >= 99 & pc < 100,'ND',
                 ifelse(pc >= 0 & pc < 1,'Normal',
                        ifelse(pc >=1 & pc < 2,'Mild',
                               ifelse(pc >= 2 & pc < 3, 'Moderate','High'))))

# Syncytia
pc <- multiplex$Syncytia
cats35 <- ifelse(pc >= 99 & pc < 100,'ND',
                 ifelse(pc >= 0 & pc < 1,'Normal',
                        ifelse(pc >=1 & pc < 2,'Mild',
                               ifelse(pc >= 2 & pc < 3, 'Moderate','High'))))

# categorizar amostras pelo CD3 pattern
tcell <- multiplex$CD3_Infiltration
cats6 <- ifelse(tcell >= 0 & tcell < 1, 'ND', 
                ifelse(tcell >= 1 & tcell < 2, 'Low infiltration',
                       ifelse(tcell >= 2 & tcell < 3, 'Moderate infiltration', 'High infiltration')))

# categorizar amostras pelo CD20 pattern
bcell <- multiplex$CD20_Infiltration
cats7 <- ifelse(bcell >= 0 & bcell < 1, 'ND',
                ifelse(bcell >= 1 & bcell < 2, 'Normal',
                       ifelse(bcell >= 2 & bcell < 3, 'Low infiltration', 
                              ifelse(bcell >= 3 & bcell < 4, 'Moderate infiltration', 'High infiltration'))))

# categorizar amostras pelo CD68 pattern
mac <- multiplex$CD68_Infiltration
cats5 <- ifelse(mac >= 0 & mac < 1, 'ND',
                ifelse(mac >= 1 & mac < 2, 'Low infiltration',
                ifelse(mac >= 2 & mac < 3, 'Moderate infiltration', 'High infiltration')))

# categorizar amostras pelo macrophage bronchial pattern
bronc <- multiplex$MAC387_Bronchial_Pattern
cats8 <- ifelse(bronc >= 0 & bronc < 1, 'ND',
               ifelse(bronc >= 1 & bronc < 2, 'Absent','Present'))

# categorizar amostras pelo macrophage alveolar pattern
alv <- multiplex$MAC387_Alveolar_Pattern
cats9 <- ifelse(alv >= 0 & alv < 1, 'ND',
                ifelse(alv >= 1 & alv < 2, 'Absent','Present'))

# categorizar amostras pelo macrophage interstitial pattern
int <- multiplex$MAC387_Interstitial_Pattern
cats10 <- ifelse(int >= 0 & int < 1,'ND',
                 ifelse(int >= 1 & int < 2, 'Absent','Present'))

# Collagen_Abundant_Macrophages
cam <- multiplex$Collagen_Abundant_Macrophages
cats23 <- ifelse(cam >= 99 & cam < 100,'ND',
                 ifelse(cam >= 0 & cam < 1,'Absent',
                        ifelse(cam >=1 & cam < 2,'Low',
                               ifelse(cam >= 2 & cam < 3, 'Moderate','High'))))

# fibrosis score
fibrosis <- multiplex$Fibrosis_score
cats13 <- ifelse(fibrosis >= 0 & fibrosis < 1,'ND',
                 ifelse(fibrosis > 1 & fibrosis <= 5,'Low',
                 ifelse(fibrosis >= 6 & fibrosis < 7, 'Moderate','High')))

# fibrin deposition
fibrin <- multiplex$Fibrin
cats14 <- ifelse(fibrin >= 0 & fibrin < 1,'ND',
                ifelse(fibrin >=1 & fibrin < 2,'Low',
                 ifelse(fibrin >= 2 & fibrin < 3, 'Moderate','High')))

# TypeII_Pneumocyte_Hyperplasia
at2 <- multiplex$TypeII_Pneumocyte_Hyperplasia
cats15 <- ifelse(at2 >= 0 & at2 < 1,'ND',
                 ifelse(at2 >=1 & at2 < 2,'Low',
                        ifelse(at2 >= 2 & at2 < 3, 'Moderate','High')))

# Alveolar_Septae_Thickening
at <- multiplex$Alveolar_Septae_Thickening
cats16 <- ifelse(at>= 0 & at < 1,'ND',
                 ifelse(at >=1 & at < 2,'Low',
                        ifelse(at >= 2 & at < 3, 'Moderate','High')))

# Alveolar_Emphysema
ae <- multiplex$Alveolar_Emphysema
cats17 <- ifelse(ae >= 99 & ae < 100,'ND',
                 ifelse(ae >= 0 & ae < 1,'Absent',
                 ifelse(ae >=1 & ae < 2,'Mild',
                        ifelse(ae >= 2 & ae < 3, 'Moderate','Severe'))))

# Alveolar_Oedema 
aed <- multiplex$Alveolar_Oedema 
cats18 <- ifelse(aed >= 99 & aed < 100,'ND',
                 ifelse(aed >= 0 & aed < 1,'Absent',
                        ifelse(aed >=1 & aed < 2,'Mild',
                               ifelse(aed >= 2 & aed < 3, 'Moderate','Severe'))))

# Broncho-pneumonia
bp <- multiplex$Broncho_pneumonia
cats19 <- ifelse(bp >= 99 & bp < 100,'ND',
                 ifelse(bp >= 0 & bp < 1,'Absent',
                        ifelse(bp >=1 & bp < 2,'Mild',
                               ifelse(bp >= 2 & bp < 3, 'Moderate','Severe'))))
                
# Bronchial_mucosal_oedema
bmo <- multiplex$Bronchial_Mucosal_Oedema
cats20 <- ifelse(bmo >= 99 & bmo < 100,'ND',
                 ifelse(bmo >= 0 & bmo < 1,'Absent',
                        ifelse(bmo >=1 & bmo < 2,'Mild',
                               ifelse(bmo >= 2 & bmo < 3, 'Moderate','Severe'))))

# Interstitial_oedema
io <- multiplex$Interstitial_Oedema
cats21 <- ifelse(io >= 99 & io < 100,'ND',
                 ifelse(io >= 0 & io < 1,'Absent',
                        ifelse(io >=1 & io < 2,'Mild',
                               ifelse(io >= 2 & io < 3, 'Moderate','Severe'))))

# Haemorrhages
ha <- multiplex$Haemorrhages
cats22 <- ifelse(ha >= 0 & ha < 1,'ND',
                        ifelse(ha >=1 & ha < 2,'Mild',
                               ifelse(ha >= 2 & ha < 3, 'Moderate','Severe')))


#Microtrombi
gt <- multiplex$Microtrombi
cats24 <- ifelse(gt >= 99 & gt < 100,'ND',
                 ifelse(gt >= 0 & gt < 1,'Absent',
                        ifelse(gt >=1 & gt < 2,'Low',
                               ifelse(gt >= 2 & gt < 3, 'Moderate','High'))))

#Arterial_Thrombi
AT <- multiplex$Arterial_Thrombi
cats25 <- ifelse(AT >= 99 & AT < 100,'ND',
                 ifelse(AT >= 0 & AT < 1,'Absent',
                        ifelse(AT >=1 & AT < 2,'Low',
                               ifelse(AT >= 2 & AT < 3, 'Moderate','High'))))

#Venous_Thrombi
VT <- multiplex$Venous_Thrombi
cats26 <- ifelse(VT >= 99 & VT < 100,'ND',
                 ifelse(VT >= 0 & VT < 1,'Absent',
                        ifelse(VT >=1 & VT < 2,'Low',
                               ifelse(VT >= 2 & VT < 3, 'Moderate','High'))))

#Angiogenesis_Intussusception
AI <- multiplex$Angiogenesis_Intussusception
cats27 <- ifelse(AI >= 99 & AI < 100,'ND',
                 ifelse(AI >= 0 & AI < 1,'Absent',
                        ifelse(AI >=1 & AI < 2,'Low',
                               ifelse(AI >= 2 & AI < 3, 'Moderate','High'))))

#Angiogenesis_Sprouting
AS <- multiplex$Angiogenesis_Sprouting
cats28 <- ifelse(AS >= 99 & AS < 100,'ND',
                 ifelse(AS >= 0 & AS < 1,'Absent',
                        ifelse(AS >=1 & AS < 2,'Low',
                               ifelse(AS >= 2 & AS < 3, 'Moderate','High'))))
#Vascular_Congestion
VC <- multiplex$Vascular_Congestion
cats29 <- ifelse(VC >= 99 & VC < 100,'ND',
                 ifelse(VC >= 0 & VC < 1,'Absent',
                        ifelse(VC >=1 & VC < 2,'Low',
                               ifelse(VC >= 2 & VC < 3, 'Moderate','High'))))

# categorizar amostras pelo length
length <- multiplex$Length
cats11 <- ifelse(length <= 1, 'Early death', 'Late death')

# categorizar by day of sampling
day <- multiplex$Day
cats12 <- ifelse(day > 0 & day <= 1, 'Day 1',
                 ifelse(day > 1 & day <= 3, 'Day 3', 
                        ifelse(day > 3 & day <= 5, 'Day 5',
                               ifelse(day > 5 & day <= 7, 'Day 7', 
                                      ifelse(day > 7 & day <= 11, 'Day 11','Day 14')))))


# removing columns for the heatmap
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("BlockID","Age", "BMI", "Weight", "Hemoglobin",	"Hematocrit", "Lymphocytes", 	"Neutrophils", 	"NCLR",	"Platelets",	"ALT",	"Glucose",	"Creatinine",	"Urea",	"LDH", "CRP",	"Temperature",
                                                "HR",	"RR",	"SpO2", "DOSAA",	"DOSAPA",	"DDOUD","Score", "Histopattern_lung", "DAD_phase","IHC_Virus_Lung", "ISH_Virus_Lung", "Cellularity", "Monocytes_Macrophages_HE", "Neutrophils_HE", 
                                               "Lymphocytes_HE",  "Plasma_cells", "Syncytia", "CD68_Infiltration", "CD3_Infiltration", "CD20_Infiltration", "MAC387_Bronchial_Pattern", "MAC387_Alveolar_Pattern", "MAC387_Interstitial_Pattern", "Collagen_Abundant_Macrophages",
                                               "Fibrosis_score", "Fibrin", "TypeII_Pneumocyte_Hyperplasia", "Alveolar_Septae_Thickening", "Alveolar_Emphysema", "Alveolar_Oedema", "Broncho_pneumonia", "Bronchial_Mucosal_Oedema", "Interstitial_Oedema", "Haemorrhages", "Microtrombi",
                                               "Arterial_Thrombi", "Venous_Thrombi", "Angiogenesis_Intussusception", "Angiogenesis_Sprouting", "Vascular_Congestion", "Length", "Category", "Day", "Day_PCA", "DOS_PCA","h.e_pattern", "DAD","Spike_IHC_pattern","ISH_pattern", "CD68_pattern", "CD3_pattern", "CD20_pattern",
                                               "Bronchial_pattern", "Alveolar_pattern", "Interstitial_pattern", "Macrophage_cluster","Score", "Patient")]

# removing columns for the correlation analysis
select_cols <- colnames(multiplex)
select_cols <- select_cols[!select_cols %in% c("BlockID", "DOSAPA", "Length", "Day", "DAD", "Category", "h.e_pattern", "Spike_IHC_pattern","ISH_pattern", "CD68_pattern", "CD3_pattern", "CD20_pattern", 
                                               "Bronchial_pattern", "Alveolar_pattern", "Interstitial_pattern", "Macrophage_cluster", "DAD.1", "Megakaryocytes_HE", "Patient")]

multiplex <- multiplex[,select_cols]
# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID
# Prepare data
mydata <- scale(df_to_cluster)

# Fazendo o heatmap
library(ComplexHeatmap)
library(circlize)
## display a palettes simultanoeusly
#options of colors for heatmaps
library(RColorBrewer)
display.brewer.all(n=10, exact.n=FALSE)
col = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(3)
col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), brewer.pal(n = 5, name = "YlOrRd"), space = "RGB")
colorRamp2(c(-0.5, 0, 0.5), c("blue", "yellow", "red"), space = "RGB"),

#heatmap virus patterns/macrophage patterns/lymphocyte patterns

#Histopattern_lung
rannot <- rowAnnotation(df = data.frame(cats), 
                        col = list(cats = c("ND" = "grey", "Acute_severe_fibrin_high_neutros" = "#a65c85ff", "Subacute_mononuclear_inflammatory_cells" = "salmon", "Haemorrhagic" = "#042333ff")), 
                        annotation_width = unit(c(0.5, 0.5), "cm"))
#DAD phase
rannot2 <- rowAnnotation(df = data.frame(cats2), 
                          col = list(cats2 = c("ND" = "grey", "Exudative" = "salmon", "Proliferative" = "#a65c85ff", "Fibrotic" = "black")), 
                          annotation_width = unit(c(0.5, 0.5), "cm"))
#spike IHC lung
rannot3 <- rowAnnotation(df = data.frame(cats3), 
                         col = list(cats3 = c("ND" = "grey" ,"Negative" = "#a65c85ff", "Positive" = "black")),
                         annotation_width = unit(c(0.5, 0.5), "cm"))
#virus ISH lung
rannot4 <- rowAnnotation(df = data.frame(cats4), 
                         col = list(cats4 = c("ND" = "grey", "Negative" = "#a65c85ff", "Positive" = "black")),
                         annotation_width = unit(c(0.5, 0.5), "cm"))

#Cellularity
rannot30 <- rowAnnotation(df = data.frame(cats30), 
                          col = list(cats30 = c("ND" = "grey", "Normal" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#Monocytes_Macrophages_HE
rannot31 <- rowAnnotation(df = data.frame(cats31), 
                          col = list(cats31 = c("ND" = "grey", "Normal" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#Neutrophils_HE
rannot32 <- rowAnnotation(df = data.frame(cats32), 
                          col = list(cats32 = c("ND" = "grey", "Normal" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Lymphocytes_HE
rannot33 <- rowAnnotation(df = data.frame(cats33), 
                          col = list(cats33 = c("ND" = "grey", "Normal" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Plasma cells
rannot34 <- rowAnnotation(df = data.frame(cats34), 
                          col = list(cats34 = c("ND" = "grey", "Normal" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Syncytia
rannot35 <- rowAnnotation(df = data.frame(cats35), 
                          col = list(cats35 = c("ND" = "grey", "Normal" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#CD3 pattern
rannot6 <- rowAnnotation(df = data.frame(cats6), 
                         col = list(cats6 = c("ND" = "grey","Low infiltration" = "salmon", "Moderate infiltration" = "#a65c85ff", "High infiltration" = "black")),
                         annotation_width = unit (c(0.5, 0.5), "cm"))
#CD20 pattern
rannot7 <- rowAnnotation(df = data.frame(cats7), 
                         col = list(cats7 = c("ND" = "grey","Normal" = "#efe350ff", "Low infiltration" = "salmon", "Moderate infiltration" = "#a65c85ff", "High infiltration" = "black")),
                         annotation_width = unit(c(0.5, 0.5), "cm"))

#CD68 pattern
rannot5 <- rowAnnotation(df = data.frame(cats5), 
                         col = list(cats5 = c("ND" = "grey", "Low infiltration" = "#efe350ff", "Moderate infiltration" = "#a65c85ff", "High infiltration" = "black")),
                         annotation_width = unit(c(0.5, 0.5), "cm"))

#Macrophage Bronchial pattern
rannot8 <- rowAnnotation(df = data.frame(cats8), 
                         col = list(cats8 = c("ND" = "grey", "Absent" = "salmon", "Present" = "black")),
                         annotation_width = unit(c(0.5, 0.5), "cm"))

#Macrophage Alveolar pattern
rannot9 <- rowAnnotation(df = data.frame(cats9), 
                         col = list(cats9 = c("ND" = "grey", "Absent" = "salmon", "Present" = "black")),
                         annotation_width = unit(c(0.5, 0.5), "cm"))

#Macrophage Interstitial pattern
rannot10 <- rowAnnotation(df = data.frame(cats10), 
                         col = list(cats10 = c("ND" = "grey", "Absent" = "salmon", "Present" = "black")),
                         annotation_width = unit(c(0.5, 0.5), "cm"))

# Collagen_Abundant_Macrophages
rannot23 <- rowAnnotation(df = data.frame(cats23), 
                          col = list(cats23 = c("ND" = "grey", "Absent" = "#efe350ff", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#Fibrosis score
rannot13 <- rowAnnotation(df = data.frame(cats13), 
                          col = list(cats13 = c("ND" = "grey", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#Fibrin deposition
rannot14 <- rowAnnotation(df = data.frame(cats14), 
                          col = list(cats14 = c("ND" = "grey", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# T2P hyperplasia
rannot15 <- rowAnnotation(df = data.frame(cats15), 
                          col = list(cats15 = c("ND" = "grey", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Alveolar thickening
rannot16 <- rowAnnotation(df = data.frame(cats16), 
                          col = list(cats16 = c("ND" = "grey", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Alveolar_Emphysema
rannot17 <- rowAnnotation(df = data.frame(cats17), 
                          col = list(cats17 = c("ND" = "grey", "Absent" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "Severe" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Alveolar oedema
rannot18 <- rowAnnotation(df = data.frame(cats18), 
                          col = list(cats18 = c("ND" = "grey", "Absent" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "Severe" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Broncho-pneumonia
rannot19 <- rowAnnotation(df = data.frame(cats19), 
                          col = list(cats19 = c("ND" = "grey", "Absent" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "Severe" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# ronchial_mucosal_oedema
rannot20 <- rowAnnotation(df = data.frame(cats20), 
                          col = list(cats20 = c("ND" = "grey", "Absent" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "Severe" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Interstitial_oedema
rannot21 <- rowAnnotation(df = data.frame(cats21), 
                          col = list(cats21 = c("ND" = "grey", "Absent" = "#efe350ff", "Mild" = "salmon", "Moderate" = "#a65c85ff", "Severe" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

# Haemorrhages
rannot22 <- rowAnnotation(df = data.frame(cats22), 
                          col = list(cats22 = c("ND" = "grey", "Mild" = "salmon", "Moderate" = "#a65c85ff", "Severe" = "#042333ff")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))
#Microtrombi
rannot24 <- rowAnnotation(df = data.frame(cats24), 
                          col = list(cats24 = c("ND" = "grey", "Absent" = "#efe350ff", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

##Arterial_Thrombi
rannot25 <- rowAnnotation(df = data.frame(cats25), 
                          col = list(cats25 = c("ND" = "grey", "Absent" = "#efe350ff", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#Venous_Thrombi
rannot26 <- rowAnnotation(df = data.frame(cats26), 
                          col = list(cats26 = c("ND" = "grey", "Absent" = "#efe350ff", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#Angiogenesis_Intussusception
rannot27 <- rowAnnotation(df = data.frame(cats27), 
                          col = list(cats27 = c("ND" = "grey", "Absent" = "#efe350ff", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#Angiogenesis_Sprouting
rannot28 <- rowAnnotation(df = data.frame(cats28), 
                          col = list(cats28 = c("ND" = "grey", "Absent" = "#efe350ff", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))
#Vascular_Congestion
rannot29 <- rowAnnotation(df = data.frame(cats29), 
                          col = list(cats29 = c("ND" = "grey", "Absent" = "#efe350ff", "Low" = "salmon", "Moderate" = "#a65c85ff", "High" = "black")),
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#length
rannot11 <- rowAnnotation(df = data.frame(cats11), 
                          col = list(cats11 = c("Early death" = "salmon", "Late death" = "light blue")), 
                          annotation_width = unit(c(0.5, 0.5), "cm"))

#day of sampling
rannot12 <- rowAnnotation(df = data.frame(cats12), 
                          col = list(cats12 = c("Day 1" = "yellow","Day 3" = "orange", "Day 5" = "salmon", "Day 7" = "red", "Day 11" = "purple", "Day 14" = "black")), 
                          annotation_width = unit(c(0.5, 0.5), "cm"))


heat <- Heatmap(mydata, cluster_columns = TRUE, cluster_rows = FALSE,
                col = colorRamp2(c(1, 0.6, 0, -0.6, -1), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                heatmap_legend_param = list(color_bar = "continuous"),
                show_row_dend = TRUE, show_column_dend = TRUE,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8))
print(heat)

print(heat+rannot)
print(heat+rannot11+rannot12)
print(rannot+rannot2+rannot3+rannot4+rannot30+rannot31+rannot32+rannot33+rannot34)
print(rannot35+rannot6+rannot7+rannot5+rannot8+rannot9+rannot10+rannot23)
print(rannot13+rannot14+rannot15+rannot16+rannot17+rannot18+rannot19+rannot20)
print(rannot21+rannot22+rannot24+rannot25+rannot26+rannot27+rannot28+rannot29)


col = colorRamp2(c(-1,-0.5, 0, 0.5, 1), magma(5, direction=-1))
col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), viridis(256) 

cmeth <- c('ward.D', 'ward.D2', "complete", "single", "average", 
           "mcquitty", 'median', "centroid")
dmeth <- c("euclidean", "manhattan", "minkowski", "pearson")

set.seed(5) 
fit <- kmeans(mydata, 3) # 3 cluster solution
pdf('heatmap_covid_fatal_clusters.pdf')
    for(cme in cmeth){
      for(dme in dmeth){
    heat <- Heatmap(mydata, cluster_columns = TRUE, cluster_rows = TRUE, clustering_distance_columns = dmeth,
                    clustering_method_columns = cmeth, clustering_distance_rows = dmeth,
                    clustering_method_rows = cmeth,
                    name = paste0(cmeth, ' ', dmeth), km = 3, 
                    col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), brewer.pal(n = 5, name = "YlOrRd"), space = "RGB"),
                    heatmap_legend_param = list(color_bar = "continuous"),
                    show_row_dend = TRUE, show_column_dend = TRUE,
                    row_names_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 8))
    
    print(heat+rannot+rannot2+rannot3+rannot4+rannot5+rannot9+rannot10+rannot6+rannot7+rannot8+rannot11
)
  }
}
dev.off()

#to print heatmap wth the number of the patients
set.seed(6)
fit <- kmeans(mydata, 3) # 3 cluster solution
pdf('heatmap_covid6.pdf')
for(cme in cmeth){
  for(dme in dmeth){
    heat <- Heatmap(mydata, cluster_columns = TRUE, cluster_rows = TRUE, clustering_distance_columns = dme,
                    clustering_method_columns = cme, clustering_distance_rows = dme,
                    clustering_method_rows = cme,
                    name = paste0(cme, ' ', dme), km = 3, 
                    col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), brewer.pal(n = 5, name = "YlOrRd"), space = "RGB"),
                    heatmap_legend_param = list(color_bar = "continuous"),
                    show_row_dend = TRUE, show_column_dend = TRUE,
                    row_names_gp = gpar(fontsize = 8),
                    column_names_gp = gpar(fontsize = 8))
    print(heat+rannot+rannot2+rannot3+rannot4+rannot5+rannot9+rannot10+rannot6+rannot7+rannot8+rannot11+rannot12+rannot15+rannot16+rannot17)
  }
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




#plot corr above with corrplot
# plot das correlacoes de acordo com o p-value
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))
library(corrplot)
cor_matrix <- corrplot(corr, method = "square", p.mat = p.mat, 
                        type = "lower", is.corr=T, sig.level = c(0.05),
                        insig = 'blank', tl.cex=0.3, tl.col = "black", number.cex=0.25, number.font=0.25, diag=F,cl.cex = 0.5, cl.ratio = 0.1,
                        mar=c(0.05,0.05,0.05,0.05), col = col4(200), order = "hclust", hclust.method = "ward.D", addgrid.col="light grey")

library(Hmisc)
cor <- rcorr(mydata, type = c("spearman")) # rcorr Calcula o p-value das correlacoes
cor$P
# plot das correlacoes de acordo com o p-value
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))
library(corrplot)
cor_matrix2 <- corrplot(cor$r, method = "square", p.mat = cor$P, 
                        type = "full", is.corr=T, insig = 'blank', sig.level = c(0.05),
                        tl.cex=0.5, tl.col = "black", number.cex=0.25, number.font=0.25, diag=F,cl.cex = 0.5, cl.ratio = 0.1,
                        mar=c(0.05,0.05,0.05,0.05), col = col4(200), addgrid.col="light grey")

#order = "hclust", hclust.method = "ward.D"
#insig = 'blank',

#Extract r and p values to plot a corr plot for the correlation of specific variables
cor$r 
write.csv(cor$r , "Correlation_spearman_rcorr_clinical_luminex_histo_early_late_death_only_last_sample.csv")
cor$P
write.csv(cor$P , "Correlation_spearman_rcorr_clinical_luminex_histo_early_late_death_only_last_sample_pvalues.csv")

#plot the correlation of specific variables using the files above
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))
library(corrplot)
cor_matrix3 <- corrplot(multiplex, method = "square", p.mat = multiplex2, 
                        type = "full", is.corr=T, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.5,pch.col = 'white',
                        tl.cex=0.5, tl.col = "black", number.cex=0.25, number.font=0.25, diag=F,cl.cex = 0.5, cl.ratio = 0.1,
                        mar=c(0.05,0.05,0.05,0.05), col = col4(200), addgrid.col="light grey")

#order = "hclust", hclust.method = "ward.D"
#insig = 'blank',

multiplex3 <- t(multiplex)
multiplex4 <- t(multiplex2)

#plot the correlation of specific variables using the files above
col4 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D",  "#B2182B", "#67001F"))
library(corrplot)
cor_matrix4 <- corrplot(multiplex, method = "square", p.mat = multiplex2, 
                        type = "full", is.corr=T, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 2,pch.col = 'white',
                        tl.cex=0.75, tl.col = "black", number.cex=0.5, number.font=0.5, diag=F,cl.cex = 0.5, cl.ratio = 0.1,
                        mar=c(0.05,0.05,0.05,0.05), col = col4(200), addgrid.col="light grey")

