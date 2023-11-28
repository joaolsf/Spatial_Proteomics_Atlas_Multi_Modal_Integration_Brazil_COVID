## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library("FactoMineR")
library("factoextra")
library("missMDA")
library(RColorBrewer)
library(ggrepel)
library(viridis)
library(ggpubr)
options(ggrepel.max.overlaps = 100)

multiplex <- read.csv("Input_SDA_Complete.csv")
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

# Prepare data
mydata <- (df_to_cluster)
mydata_luminex <- mydata[,1:72]
mydata_clinical <- mydata[,73:89]
mydata_IMC <- mydata[,90:202]

# Perform the (regularized) iterative MFA algorithm with the number of dimensions selected in the previous step, using the function imputePCA:
mydata2 <- imputeMFA(mydata[,1:202], group = c(72, 17, 113), #luminex, clinical, IMC+histo continuous, histo categorical
                     type = c("s", "s", "s"),
                     ncp = 5)
mydata2 <- mydata2$completeObs
write.csv(mydata2, "Input_SDA_Complete.csv")


day1 <- read.csv("Day1.csv")
day3 <- read.csv("Day3.csv")
day5 <- read.csv("Day5.csv")
day7 <- read.csv("Day7.csv")
day11 <- read.csv("Day11.csv")
day14 <- read.csv("Day14.csv")

day1 <- as.data.frame(day1)
day1 <- day1[,-1]
day3 <- as.data.frame(day3)
day3 <- day3[,-1]
day5 <- as.data.frame(day5)
day5 <- day5[,-1]
day7 <- as.data.frame(day7)
day7 <- day7[,-1]
day11 <- as.data.frame(day11)
day11 <- day11[,-1]
day14 <- as.data.frame(day14)
day14 <- day14[,-1]

# Create array adding a 4th dimention for the number of approaches
array <- array(data = c(unlist(day1), unlist(day3), unlist(day5), unlist(day7), unlist(day11), unlist(day14)),
               dim = c(16, 202, 6, 3), dimnames = list(rownames(day1),
                                                        colnames(day1)))

# Create array for 16 samples, 202 variables and 6 time points
array2 <- array(data = c(unlist(day1), unlist(day3), unlist(day5), unlist(day7), unlist(day11), unlist(day14)),
               dim = c(16, 202, 6), dimnames = list(rownames(day1),
                                                       colnames(day1)))




## ----setup--------------------------------------------------------------------
# #Load the library SDA4D
# #you can use devtools, for example
#if(!require(devtools)){install.packages('devtools')}
#devtools::install_github('https://github.com/marchinilab/SDA4D')
# #load library
library(SDA4D)

## ----generatedata-------------------------------------------------------------
#set.seed(42)
#generatedData <- generate_example_data()
#lapply(generatedData,dim)
#write.csv(generatedData$dataTensor, "dataTensor.csv")
#dim(generatedData$dataTensor)

## ----runmethod,results=FALSE--------------------------------------------------
sda <- RunSDA4D(data_tensor = array,
               dimn_vector = c(16,202,6,3),
               max_iters = 2000,
               num_components = 4,
               stopping=FALSE)

saveRDS(sda, "SDA_model.rds")
sda <- readRDS("SDA_model.rds")
## ----plotELBO,fig.width=8,fig.height=6,echo=FALSE-----------------------------
#res<-list(Neg_FE=exp(rnorm(100)))
  #plot(-log(-res$Neg_FE[-1]),col='red',pch=4,xlab='Iteration',main='log(ELBO) over 100 iterations')
  #lines(-log(-res$Neg_FE[-1]),col='red')

## ----plotELBOqplot,echo=FALSE,fig.width=6,fig.height=4,message = FALSE,warning = FALSE----
suppressWarnings(library(ggplot2))
startfrom=5
ELBO_plot<-qplot(x=c(startfrom:length(sda$ELBO)),
                 y=sda$ELBO[-c(1:(startfrom-1))],
                 geom=c("point", "line"))+
            ylab('ELBO')+
            xlab('Iteration')+
            ggtitle('ELBO over 100 iterations')
print(ELBO_plot)

## ----outputnames--------------------------------------------------------------
names(sda)
sda$maximumiteration
length(sda$ELBO)

## ----outputnamesanddims-------------------------------------------------------
lapply(sda$A,dim) # A - looks like the values of each component (factor) per patient
factor_loadings_patient <- as.data.frame(sda$A$mu)
write.csv(factor_loadings_patient, "factor_loadings_patient.csv")

lapply(sda$B,dim) # B- scores of each approach per component (factor)
factor_loadings_approach <- as.data.frame(sda$B$nu)
write.csv(factor_loadings_approach, "factor_loadings_approach.csv")

lapply(sda$D,dim) # D- scores of each component (factor) over time
factor_loadings_time <- as.data.frame(sda$D$alpha)
write.csv(factor_loadings_time, "factor_loadings_time.csv")

lapply(sda$WS,dim) #WS - looks like the values of each component per parameter
factor_loadings_parameters <- as.data.frame(sda$WS$mom1)
write.csv(factor_loadings_parameters, "factor_loadings_parameters.csv")

factor_loadings_parameters <- as.data.frame(sda$WS$m)
write.csv(factor_loadings_parameters, "factor_loadings_parameters_3.csv")

identical(sda$WS$mom1,sda$WS$m*sda$WS$gamma)

library(SDAtools)

# Simulate data
set.seed(42)
data <- simulate_2D_data()
export_data(data$Y, name = "simulated.data", path = "./data-raw/")
write.csv(data$Y, "simulate_data.csv")

run_SDA(out = "./data-raw/simulation_results2",
        data = "./data-raw/simulated.data",
        max_iter = 200,
        save_freq = 200)


















