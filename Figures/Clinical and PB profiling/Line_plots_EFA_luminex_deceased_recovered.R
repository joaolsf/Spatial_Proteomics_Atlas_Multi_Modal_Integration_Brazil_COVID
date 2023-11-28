library(ggplot2)
library(reshape)
library(reshape2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(tidyr)
library(forcats)
library(Rmisc)
library(ggpubr)


# reading file (importar o arquivo)
multiplex <- read.csv("COVID_Luminex_deceased_recovered_factor_scores.csv")

longer_data <- multiplex %>%
  pivot_longer(Factor_1:Factor_5, names_to = "question", values_to = "response")

multiplex %>%
  pivot_longer(Factor_1:Factor_5, names_to = "question", values_to = "response") %>%
  ggplot(aes(x =Day_of_Sampling, y=response, group=Category)) + 
  geom_errorbar(aes(ymin=Factor_1:Factor_5-sd, ymax=Factor_1:Factor_5+sd), width=.1) +
  geom_line(aes(color=Category), size = 0.5) + geom_point(aes(color=Category)) +
  facet_wrap(~question, scale="free", ncol = 5) +
  labs(x = "Day of Sampling") + scale_color_manual(values=c("salmon", "light blue", "dark grey")) + theme_bw(base_size = 10)


"#6b4596ff", "#de7065ff", "#f7cb44ff"
"#440154FF", "#29AF7FFF", "#f7cb44ff"

#line plot for each signature

#The function below will be used to calculate the mean and the standard deviation, for the variable of interest, in each group :
  
  #+++++++++++++++++++++++++
  # Function to calculate the mean and the standard deviation
  # for each group
  #+++++++++++++++++++++++++
  # data : a data frame
  # varname : the name of a column containing the variable
  #to be summariezed
  # groupnames : vector of column names to be used as
  # grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#Summarize your dataset
#Standard deviation calculation for each factor 
df1 <- data_summary(multiplex, varname="Factor_1", 
                    groupnames=c("Category", "Day_of_Sampling"))
head(df1)

df2 <- data_summary(multiplex, varname="Factor_2", 
                    groupnames=c("Category", "Day_of_Sampling"))
head(df2)

df3 <- data_summary(multiplex, varname="Factor_3", 
                    groupnames=c("Category", "Day_of_Sampling"))
head(df3)

df4 <- data_summary(multiplex, varname="Factor_4", 
                    groupnames=c("Category", "Day_of_Sampling"))
head(df4)

df5 <- data_summary(multiplex, varname="Factor_5", 
                    groupnames=c("Category", "Day_of_Sampling"))
head(df5)

#Standard error or mean calculation for each factor 
df1_1 <- summarySE(multiplex, measurevar="Factor_1", groupvars=c("Category", "Average_DOS"))
head(df1_1)

df2_1 <- summarySE(multiplex, measurevar="Factor_2", groupvars=c("Category", "Average_DOS"))
head(df2_1)

df3_1 <- summarySE(multiplex, measurevar="Factor_3", groupvars=c("Category", "Average_DOS"))
head(df3_1)

df4_1 <- summarySE(multiplex, measurevar="Factor_4", groupvars=c("Category", "Average_DOS"))
head(df4_1)

df5_1 <- summarySE(multiplex, measurevar="Factor_5", groupvars=c("Category", "Average_DOS"))
head(df5_1)

write.csv(df1_1, "Factor1_means.csv")
write.csv(df2_1, "Factor2_means.csv")
write.csv(df3_1, "Factor3_means.csv")
write.csv(df4_1, "Factor4_means.csv")
write.csv(df5_1, "Factor5_means.csv")

#Plot line plots with SEM values
compare_means(Factor_1 ~ Category, data = multiplex, 
              group.by = "Day_of_Sampling")


df1_1$title <- "PB signature 1"
ggplot(df1_1, aes (x=Average_DOS, y=Factor_1, group=Category)) +
  geom_errorbar(aes(ymin=Factor_1-se, ymax=Factor_1+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Symptoms") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 1") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df2_1$title <- "PB signature 2"
ggplot(df2_1, aes (x=Average_DOS, y=Factor_2, group=Category)) +
  geom_errorbar(aes(ymin=Factor_2-se, ymax=Factor_2+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Symptoms") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 2") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df3_1$title <- "PB signature 3"
ggplot(df3_1, aes (x=Average_DOS, y=Factor_3, group=Category)) +
  geom_errorbar(aes(ymin=Factor_3-se, ymax=Factor_3+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Symptoms") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 3") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df4_1$title <- "PB signature 4"
ggplot(df4_1, aes (x=Average_DOS, y=Factor_4, group=Category)) +
  geom_errorbar(aes(ymin=Factor_4-se, ymax=Factor_4+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Symptoms") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 4") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df5_1$title <- "PB signature 5"
ggplot(df5_1, aes (x=Average_DOS, y=Factor_5, group=Category)) +
  geom_errorbar(aes(ymin=Factor_5-se, ymax=Factor_5+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Symptoms") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 5") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

#Standard error or mean calculation for each factor 
df1_1 <- summarySE(multiplex, measurevar="Factor_1", groupvars=c("Category", "Day_of_Sampling"))
head(df1_1)

df2_1 <- summarySE(multiplex, measurevar="Factor_2", groupvars=c("Category", "Day_of_Sampling"))
head(df2_1)

df3_1 <- summarySE(multiplex, measurevar="Factor_3", groupvars=c("Category", "Day_of_Sampling"))
head(df3_1)

df4_1 <- summarySE(multiplex, measurevar="Factor_4", groupvars=c("Category", "Day_of_Sampling"))
head(df1_1)

df5_1 <- summarySE(multiplex, measurevar="Factor_5", groupvars=c("Category", "Day_of_Sampling"))
head(df2_1)

df1_1$title <- "PB signature 1"
ggplot(df1_1, aes (x=Day_of_Sampling, y=Factor_1, group=Category)) +
  geom_errorbar(aes(ymin=Factor_1-se, ymax=Factor_1+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Hospitalisation") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 1") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df2_1$title <- "PB signature 2"
ggplot(df2_1, aes (x=Day_of_Sampling, y=Factor_2, group=Category)) +
  geom_errorbar(aes(ymin=Factor_2-se, ymax=Factor_2+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Hospitalisation") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 2") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df3_1$title <- "PB signature 3"
ggplot(df3_1, aes (x=Day_of_Sampling, y=Factor_3, group=Category)) +
  geom_errorbar(aes(ymin=Factor_3-se, ymax=Factor_3+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Hospitalisation") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 3") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df4_1$title <- "PB signature 4"
ggplot(df4_1, aes (x=Day_of_Sampling, y=Factor_4, group=Category)) +
  geom_errorbar(aes(ymin=Factor_4-se, ymax=Factor_4+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Hospitalisation") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 4") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)

df5_1$title <- "PB signature 5"
ggplot(df5_1, aes (x=Day_of_Sampling, y=Factor_5, group=Category)) +
  geom_errorbar(aes(ymin=Factor_5-se, ymax=Factor_5+se, color=Category), size = 0.75, width=0.75, position=position_dodge(0.05)) +
  geom_line(aes(color=Category), size = 2) + geom_point(aes(color=Category), size=10, shape=20) +
  xlab("Days of Hospitalisation") + 
  ylab("Sample loading score") +
  ggtitle("PB signature 5") +
  expand_limits(y=0) +                        # Expand y range
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 18, base_line_size = 0.75) +
  theme(legend.justification=c(),
        legend.position='bottom',
        axis.text.y = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=18, angle=0, vjust=1, hjust=0.5),
        strip.text.x = element_text(size = 24, colour = "white", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#800080", size=1.5, linetype="blank")) +              # Position legend in bottom right
  scale_color_manual(values=c("salmon", "light blue", "dark grey")) + facet_grid(. ~ title)


