library(reshape)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(tidyr)
library(forcats)
library(Rmisc)


#BOX PLOTS METHOD clinical parameter
multiplex <- read.csv("COVID_Clinical_deceased_recovered_boxplots.csv")

multiplex$title <- "Days in Mechanical Ventilation"
p <- ggplot(multiplex, aes(x=Category, y=DMV, fill=Category)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.shape = 2)  + 
  geom_jitter(size = 4, shape=16, position=position_jitter(0.1)) +
  labs(y = "Days") +
  theme_bw(base_size = 24) +
  theme(legend.justification=c(),
        legend.position='bottom', plot.title = element_text(size=20), 
        axis.text.y = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5), 
        strip.text.x = element_text(size = 20, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="black", size=1.5, linetype="blank"))

p + scale_fill_manual(values = c("dark blue", "light grey")) + scale_x_discrete(limits=c("Recovered", "Fatal"))  + facet_grid(. ~ title) #+ scale_y_continuous(limit = c(0, 40))

#Box plot per each PB biomarker
multiplex <- read.csv("COVID_Clinical_deceased_recovered_hierarchical_clusters_boxplots.csv")

multiplex$title <- "Syndecan-1"
p <- ggplot(multiplex, aes(x=Category, y=Syndecan1, fill=Category)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.shape = 2)  + 
  geom_jitter(size = 4, shape=16, position=position_jitter(0.1)) +
  labs(y = "Syndecan-1") +
  theme_bw(base_size = 24) +
  theme(legend.justification=c(),
        legend.position='bottom', plot.title = element_text(size=20), 
        axis.text.y = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5), 
        strip.text.x = element_text(size = 20, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank"))

p + scale_fill_manual(values = c("light blue", "salmon")) + scale_x_discrete(limits=c("Cluster 1", "Cluster 2"))  + facet_grid(. ~ title) #+ scale_y_continuous(limit = c(0, 40))

multiplex$title <- "Days of Disease Onset until Death"
p <- ggplot(multiplex, aes(x=Category, y=DDOUD, fill=Category)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.shape = 2)  + 
  geom_jitter(size = 4, shape=16, position=position_jitter(0.1)) +
  labs(y = "Days") +
  theme_bw(base_size = 24) +
  theme(legend.justification=c(),
        legend.position='bottom', plot.title = element_text(size=20), 
        axis.text.y = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5), 
        strip.text.x = element_text(size = 18, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="black", size=1.5, linetype="blank"))

p + scale_fill_manual(values = c("light blue", "salmon")) + scale_x_discrete(limits=c("Cluster 1", "Cluster 2"))  + facet_grid(. ~ title) #+ scale_y_continuous(limit = c(0, 40))

# #800080" = purple
# #440154FF" = deep purple
# #F79C79 = orange
# #21908CFF" = green
# #0080FF - blue

#BAR PLOT
multiplex$title <- "Days of Disease Onset until Death"
p2 <- ggbarplot(multiplex, x='Category', y='DDOUD', fill='Category', width = 0.5, add = c("median_iqr"), 
                label = FALSE, error.plot = "linerange") +
  geom_jitter(size = 4, shape=16, position=position_jitter(0.1)) +
  labs(title="DDOUD", y = "Days") +
  theme_bw(base_size = 24) +
  theme(legend.justification=c(),
        legend.position='bottom', plot.title = element_text(size=20), 
        axis.text.y = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=24, angle=0, vjust=1, hjust=0.5), 
        strip.text.x = element_text(size = 18, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="black", size=1.5, linetype="blank"))

p2 + scale_fill_manual(values = c("light blue", "salmon")) + scale_x_discrete(limits=c("Cluster 1", "Cluster 2")) + facet_grid(. ~ title) #+ scale_y_continuous(limit = c(0, 35))

# #EFC000FF - yellow
# #0073C2FF - blue

#Box plot kmeans clusters
multiplex <- read.csv("COVID_clinical_data_kmeans_clusters.csv")

multiplex$title <- "Days of Disease Onset until Death"
p3 <- ggbarplot(multiplex, x='kmeans_cluster', y='DDOUD', fill='kmeans_cluster', width = 0.5, add = c("median_iqr"), 
                label = FALSE, error.plot = "linerange") +
  geom_jitter(size=4, shape=16, position=position_jitter(0.1)) +
  labs(title="DDOUD",y = "Days") +
  theme_bw(base_size = 22) +
  theme(legend.justification=c(),
        legend.position='bottom', plot.title = element_text(size=20), 
        axis.text.y = element_text(face = "plain", size=22, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=16, angle=0, vjust=1, hjust=0.5), 
        strip.text.x = element_text(size = 18, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="black", size=1.5, linetype="blank"))

p3 + scale_fill_manual(values = c("#EFC000FF","#0073C2FF")) + scale_x_discrete(limits=rev(c("Kmeans cluster 1", "Kmeans cluster 2"))) + facet_grid(. ~ title) 

scale_x_discrete(limits = ))
