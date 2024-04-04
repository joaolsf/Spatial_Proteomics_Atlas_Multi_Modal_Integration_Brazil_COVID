
# 1- Load libraries
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(tidyverse)
library(cowplot)
library(reshape2)
library(MOFA2)

# This option is to be used with use_basilisk = FALSE)
library(reticulate)
#py_install("mofapy2", envname = "mofa_env2", method="conda", pip=TRUE) #just run if you don't have this yet
use_condaenv(condaenv="mofa_env2", required=TRUE)
mofa2 <- import("mofapy2")
# py_list_attributes(mofa2) #this is important to check whether the module contains all required attributes

# 2- Preprocessing

# 2.1- Input missing data with missMDA package (PCA approach)
library(missMDA)
library(FactoMineR)
library(factoextra)

multiplex <- read.csv("COVID_early_late_death_MOFA_MEFISTO_missing_data.csv")

# data frame para cluster
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID

#subsetting active variables for the PCA
df_to_cluster.active <- df_to_cluster[,1:207] 

# Prepare data
mydata <- (df_to_cluster.active)

#Perform the (regularized) iterative PCA algorithm with the number of dimensions selected in the previous step, using the function imputePCA:
res.imp <- imputePCA(mydata, ncp = 6)

mydata_complete <- res.imp$completeObs
write.csv(mydata_complete, "COVID_early_late_death_MOFA_MEFISTO_missing_data_complete_factor_values.csv")

# 2.2- Scale data

multiplex <- read.csv("COVID_early_late_death_MOFA_MEFISTO_missing_data_complete.csv")
df_to_cluster <- multiplex[,-1]
rownames(df_to_cluster) <- multiplex$RecordID
df_to_cluster_scale <- df_to_cluster[,1:202]
# df_to_cluster_supp <-df_to_cluster[,203:207]

# Scale data
mydata <- scale(df_to_cluster_scale)
mydata2 <- cbind.data.frame(mydata, df_to_cluster_supp)
write.csv(mydata, "COVID_early_late_death_MOFA_MEFISTO_missing_data_complete_scaled.csv")

# 2.3- Prepare MOFA?MEFISTO input df
mydata3 <- read.csv("COVID_early_late_death_MOFA_MEFISTO_missing_data.csv")
library(reshape)
data.m <- melt(mydata3)
write.csv(data.m, "COVID_early_late_death_MOFA_MEFISTO_iteration8_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################

# 3- Create the MOFA object and train the model
#To use the multi-group function, this has to be defined with create_mofa function.
mydata4 <- read.csv("COVID_early_late_death_MOFA_MEFISTO_iteration8_v2.csv")

#to remove the "group" column
select_cols <- colnames(mydata4)
select_cols <- select_cols[!select_cols %in% c("group")]
mydata4 <- mydata4[,select_cols]

#Create MOFAobject
MOFAobject_untrained <- create_mofa(mydata4) 
MOFAobject_untrained
#groups = "group"

# 3.1- If using the time as covariate for MEFISTO analysis
#Next, we want to add the time information for each sample, which we can do using set_covariates. 
#As the time is already contained in the data.frame that we passed to the MOFAobject in create_mofa, 
#we can just specify to use this column as covariate. 
#Alternatively, we could also supply a new matrix or data frame providing the sample names and covariates. See ?set_covariates.
head(samples_metadata(MOFAobject_untrained))

MOFAobject_untrained <- set_covariates(MOFAobject_untrained, covariates = "Average_Day_of_Symptoms")
MOFAobject_untrained

#"Average_Day_of_Symptoms"

# Visualise data structure overview
#Visualise the number of views (rows) and the number of groups (columns) exist, 
#what are their corresponding dimensionalities 
#and how many missing information they have (grey bars).
plot_data_overview(MOFAobject_untrained)

#Show data and time covariate Before moving towards the training, 
#let’s first take a look at the data in the untrained object. 
#The right plots show the time values for each sample.
gg_input <- plot_data_overview(MOFAobject_untrained,
                               show_covariate = TRUE,
                               show_dimensions = TRUE) 
gg_input


# 3.2- Define MOFA options

#3.2.1- Data options
#Important arguments:
#scale_groups: scale groups to the same total variance? Default is FALSE
#scale_views: scale views to the same total variance? Default is FALSE
#views: views names
#groups: groups names
data_opts <- get_default_data_options(MOFAobject_untrained)
data_opts$center_groups <- FALSE
data_opts

#3.2.2- Model options
model_opts <- get_default_model_options(MOFAobject_untrained)
model_opts$num_factors <- 5 #play with this

#3.2.3- Train options
train_opts <- get_default_training_options(MOFAobject_untrained)
train_opts$seed <- 2020
train_opts$convergence_mode <- "slow" # use "fast" for faster training

#3.2.4- MEFISTO options (if using time or space as covariate)
mefisto_opts <- get_default_mefisto_options(MOFAobject_untrained) #only if using mefisto
#mefisto_opts$warping <- TRUE #only if aligning time points that are misaligned like in days of symptoms
#mefisto_opts$warping_ref <- "Mouse" # Here, we specify that we want to align the times across species using Mouse as a reference group. 
#mefisto_opts$new_values <- 1:14  #Furthermore, we can optionally specify time points (when using warping these correspond to times in the reference group), for which to interpolate the factors in all groups.
#mefisto_opts$new_values <- matrix(0:28, nrow =1) # set time points to interpolate factors to

#3.3- Prepare MOFA
#We then pass all these options to the object.
MOFAobject_untrained <- prepare_mofa(
  object = MOFAobject_untrained,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts,
  mefisto_options = mefisto_opts
) 

#3.4- Train the MOFA model 
#Now, we are ready to train the model.
outfile <- "COVID_early_late_death_MEFSTO_model_day_of_symptoms_iteration8_v2.hdf5"
MOFAobject <- run_mofa(MOFAobject_untrained, outfile = outfile)


#4- Downstream Analysis
# Load the model
outfile2 <- "COVID_early_late_death_model_iteration1_MEFISTO_day_of_symptoms.hdf5"
MOFAobject <- load_model(outfile2, load_interpol_Z = FALSE) #In addition to the factor values at the observed data points we can optionally load the interpolated data.

#Define colors for the visualisations
group.colors <- c(
  "Early death" = "salmon",
  "Late death" = "light blue")

#4.1- Imputation of missing values
#With the impute function all missing values are imputed based on the MOFA model. 
#The imputed data is then stored in the imputed_data slot of the MOFAobject and can be accessed via the get_imputed_data function.
MOFAobject <- impute(MOFAobject)

#4.2- Correlation between factors
#A good sanity check is to verify that the Factors are largely uncorrelated.
#In MOFA there are no orthogonality constraints such as in Principal Component Analysis, 
#but if there is a lot of correlation between Factors this suggests a poor model fit.
#Reasons? Perhaps you used too many factors or perhaps the normalisation is not adequate.
plot_factor_cor(MOFAobject)

#4.3- Plot variance decomposition
#4.3.1- Variance decomposition by Factor
#The most important insight that MOFA generates is the variance decomposition analysis.
#This plot shows the percentage of variance explained by each factor across each data modality (and group, if provided). 
#It summarises the sources of variation from a complex heterogeneous data set in a single figure.
plot_variance_explained(MOFAobject, max_r2=20)
plot_variance_explained(MOFAobject, factor=c(1, 2, 3, 4, 5), legend = T) #select specific factors to plot their variance explained

#4.3.2- Plot variance explained per factor across groups
#Quantifying the variance explained per factor across groups and views is probably the most important plot that MOFA+ generates. 
#It summarises the (latent) signal from a complex heterogeneous data set in a single figure.
#plot_variance_explained(MOFAobject, x="group", y="factor")
#We can also plot the total variance explained per group (with all factors) by adding the argument plot_total = TRUE. 
#Notably, only 10 factors are sufficient to capture between 35% and 55% of the transcriptional variance per embryo
plot_variance_explained(MOFAobject, x = "group", y = "factor", plot_total = T)[[2]]

#4.3.3- Plot variance explained for individual features
#We can also inspect the variance explained by the MOFA factors for individual features. 
#A high R2 implies that the MOFA factors captures most of the variation for this feature, 
#whereas small values means that the variation for this feature is not explained by the model (i.e. it is considered as noise):
features <- c("Syndecan1","ICAM1","VCAM1","PD_L1")
plot_variance_explained_per_feature(MOFAobject, view = "Luminex", features = features)

#4.3.4- Total variance explained per view
#A reasonable question is whether the model is providing a good fit to the data.
#For this we can plot the total variance explained (using all factors). The resulting values will depend on the nature of the data set, the number of samples, the number of factors, etc. Some general guidelines:
#Noisy data sets with strong non-linearities will result in small amounts of variance explained (<10%).
#The higher the number of samples the smaller the total variance explained
#The higher the number of factors, the higher the total variance explained.
#MOFA is a linear and sparse model. This is helpful to prevent overfitting, but it will never explain 100% of the variance, even if using a lot of Factors.
#In this data set, using only K=5 factors the model explains up to ~70% of the variation in the IMC and ~50% of the variation in the luminex data. This is quite remarkable for a linear model.
plot_variance_explained(MOFAobject, plot_total = T)[[2]]

#4.3.5- Characterisation of Factors
#There are a few systematic strategies to characterise the molecular etiology underlying the MOFA Factors and to relate them to the sample covariates:
#Association analysis between the sample metadata and the Factor values.
#Inspection of factor values.
#Inspection of the feature weights.

#4.3.5.1- Association analysis
#Let’s test the association between MOFA Factors and covariates:
correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Group","Day_of_Sampling","Average_Day_of_Symptoms"), 
                                  plot="log_pval")

#4.3.5.2- Inspection of factor values
#Each factor captures a different source of variability in the data. 
#Mathematically, each Factor is defined by a linear combination of the input features. 
#As the data is centered prior to running MOFA, each Factor ordinates cells along a one-dimensional axis that is centered at zero. 
#Samples with different signs manifest opposite phenotypes along the inferred axis of variation, with higher absolute value indicating a stronger effect.
#Note that the interpretation of MOFA factors is analogous to the interpretation of the principal components in PCA.
plot_factor(MOFAobject, factors = 1, color_by = "Factor1")

#Alternative graph
plot_factor(MOFAobject, 
            factor = 2, 
            color_by = "Group", 
            dot_size = 1,
            dodge = TRUE,
            stroke = 0.4,
            add_violin = T,
            add_boxplot = T) +
  scale_fill_manual(values=group.colors) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())
#Alternative graph
plot_factor(MOFAobject,
            factor = 2,
            color_by = "Group", 
            scale = TRUE, 
            add_violin = TRUE, color_violin = TRUE, 
            dodge = TRUE, dot_size = 1, legend = TRUE
) + scale_color_manual(values=group.colors) + scale_fill_manual(values=group.colors)


#MEFISTO ONLY ANALYSIS - Factor overview and visualization
# 4.3.6 Variance decomposition and factor correlation
# To obtain a first overview of the factors we can take a look at the variance that a factor explains in each view.
p <- plot_variance_explained(MOFAobject, plot_total = T)
p[[1]] + theme(axis.text.x = element_text(angle = 0))
p[[2]] + theme(axis.text.x = element_text(angle = 0))

# 4.3.7- Factors versus temporal covariates
# Factors versus temporal covariate. To investigate the inferred factors, 
# we can plot them against the days of sampling, days of symptoms and colour them by the metadata of the samples such as the disease group
# Factors colored by disease group and vs days of symptoms
plot_factors_vs_cov(MOFAobject, color_by = "Group", factor = 1, covariate = "Day_of_Sampling", scale = FALSE, show_variance = FALSE) +
  scale_fill_manual(values = group.colors) +
  stat_summary(aes(col = color_by), geom="line", fun = "mean", size=2) +
  scale_color_manual(values = group.colors) + geom_point(aes(color=color_by), size=2, shape=19) +
  scale_x_continuous(breaks=c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)) +
  theme_bw(base_size = 20) + labs(title="Integrated signature 1", y = "Value") +
  theme(plot.title = element_text(size=20), 
        axis.text.y = element_text(face = "plain", size=20, angle=0, vjust=1, hjust=0.5),
        axis.text.x = element_text(face = "plain", size=20, angle=0, vjust=1, hjust=0.5), 
        strip.text.x = element_text(size = 20, colour = "black", face = "plain"),
        strip.background = element_rect(color=NULL, fill="#E7B800", size=1.5, linetype="solid"))

# Factors colored by disease group and vs days of sampling
plot_factors_vs_cov(MOFAobject, color_by = "Group", factor = 3) +
  scale_fill_manual(values = group.colors) +
  stat_summary(aes(col = color_by), geom="line", fun = "mean", size=2) +
  scale_color_manual(values = group.colors) + geom_point(aes(color=color_by), size=6, shape=20) +
  scale_x_continuous(breaks=c(1, 3, 5, 7, 11, 14)) +
  theme_bw(base_size = 16)  

#4.3.8- Alternative graphs - when using temporal alignment
#Before alignment, the timing of common developmental process does not match, after alignment, we see a clearer match.
# Factors versus developmental time (before alignment)
plot_factors_vs_cov(MOFAobject, color_by = "Group", covariate = "Average_Day_of_Symptoms",  warped = FALSE) 
# Factors versus developmental time (after alignment)
plot_factors_vs_cov(MOFAobject, color_by = "Group", covariate = "Average_Day_of_Symptoms", scale = TRUE)

#4.3.9- Plot feature weights
#How do we interpret the weights?
#The weights provide a score for each feature on each factor. 
#Features with no association with the corresponding factor are expected to have values close to zero, 
#whereas features with strong association with the factor are expected to have large absolute values. 
#The sign of the weights indicates the direction of the effect: a positive weights indicates that the feature has higher levels in the cells with positive factor values, and vice-versa.

#4.3.9.1- Plot feature weights for each view
#By looking at the variance explained plot, we saw that Factor 1 captures variation in all data modalities. 
#Out of all omics, the somatic mutation data is a good place to start, as somatic mutations are very sparse, easy to interpret and any change in the DNA is likely to have downstream consequences to all other molecular layers. 
#Define a helper function to plot weights, it will reduce the amount of code in the vignette
plot_weights_fn <- function(MOFAobject, factor=1, view=1, nfeatures=10) {
  p1 <- plot_weights(MOFAobject, 
                     factors = factor, 
                     view = view,
                     nfeatures = nfeatures,
                     text_size = 4, dot_size = 2)
  
  p2 <- plot_top_weights(MOFAobject, 
                         factors = factor, 
                         view = view,
                         nfeatures = nfeatures)
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
  return(p)
}

plot_weights_fn(MOFAobject, factor=3, view="Luminex", nfeatures=10)
plot_weights_fn(MOFAobject, factor=3, view="Clinical", nfeatures=10)
plot_weights_fn(MOFAobject, factor=3, view="IMC", nfeatures=10)

# Get the weights to save to a csv file - top negative ones
weights_df_top_negative = get_weights(MOFAobject, factors = 1, as.data.frame = TRUE, views = c('Clinical'), scale = TRUE) %>%
  arrange(value) %>% head(20) %>%
  left_join(., features_metadata(MOFAobject), by = c("feature", "view")) %>%
  select(feature, value, view)

# Get the weights to save to a csv file - top positive ones 
weights_df_top_positive <- get_weights(MOFAobject, factors = 1, as.data.frame = TRUE, views = c('Clinical'), scale = TRUE) %>%
  arrange(-value) %>% head(20) %>%
  left_join(., features_metadata(MOFAobject), by = c("feature", "view")) %>%
  select(feature, value, view) 

#bar plots with the top feature weight values for each view
#How to Rotate Axis Labels in ggplot2 - https://www.statology.org/rotate-axis-labels-ggplot2/
weights <- read.csv("Factor3_IMC_features_weights.csv")
plot_ <- ggplot(weights,
                aes(x= reorder(Feature,
                               (Value), decreasing = FALSE), y = (Value))) +
  geom_bar(stat = "identity", 
           fill = "#21908CFF", width=0.7, position=position_dodge(width=1)) +
  theme_bw(base_size = 10) +
  coord_flip() +
  xlab("") + 
  ylab("Factor weight") +
  ggtitle("Tissue parameters - Integrated signature 3") +
  theme(plot.title = element_text(size=16), 
        axis.text.y = element_text(face = "bold", size=16, angle=0, vjust=.5, hjust=1),
        axis.text.x = element_text(face = "bold", size=16, angle=0, vjust=.5, hjust=1)) 
plot_

# #800080" = purple
# #440154FF" = deep purple
# #F79C79 = orange
# #21908CFF" = green

#views = c('IMC','Luminex','Clinical')


#Plot the weights:
plot_weights(MOFAobject,
             view = "IMC",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T, dot_size = 1)           # Scale weights from -1 to 1
#An alternative visualisation to the full distribution of weights is to do a line plot that displays only the top features with the corresponding weight sign on the right:
plot_top_weights(MOFAobject,
                 view = "Clinical",
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T)           # Scale weights from -1 to 1


#If you are not interested in the directional of the effect, you can take the absolute value of the weights (abs=TRUE). 
#You can also highlight some variables of interest using the argument manual:
plot_weights(MOFAobject,
             view = "RNA",
             factor = 1,
             nfeatures = 5,
             manual = list(c("Snai1","Mesp1","Phlda2"), c("Rhox5","Elf5")),
             color_manual = c("darkgreen","red"),
             scale = T,
             abs = F)

#4.3.9.2- Factor weights
#Next we have a look at the weights of the factor.
#Inspection of top genes per factor
#To obtain insights into the molecular signatures associated with each factor in view, we can inspect the genes that have high weights on the factors.
plot_top_weights(MOFAobject, view = "Luminex", factors = 2)
#or
plot_data_vs_cov(MOFAobject, factor = 1, view = "Luminex", color_by = "Group", shape_by = "Group", features = 6) +
  scale_fill_manual(values = group.colors)

#Alternative graphs
plot_data_vs_cov(MOFAobject, features = 3,
                 color_by = "Group") + 
  scale_fill_manual(values = group.colors)

#Plot the Factor values and color based on specific features.
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Group",
            add_violin = TRUE,
            dodge = TRUE)

#We can also plot Factor values coloured by other covariates, for example day of sampling
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Day_of_Sampling",
            dodge = TRUE,
            add_violin = TRUE)

#Visualization of gene expression patterns given by the Factors
#We expect negative weights For Factor 1 to be marker genes of ExE Endoderm. 
#If we plot Factor 1 coloring cells by mRNA expression of the genes with the largest weights:
genes <- c("Ttr","Apom")
for (i in genes) {
  p <- plot_factor(MOFAobject,
                   factor = 1,
                   color_by = i
  ) + scale_colour_gradientn(colours = terrain.colors(10))
  print(p)
}

#4.3.9.2- Plot molecular signatures in the input data
#In this case we have a large amount of genes that have large positive and negative weights.
#The function plot_data_scatter generates a scatterplot of Factor 1 values (x-axis) versus expression values (y-axis) for the top 4 markers with largest positive weight. 
#Samples are colored by disease group:
plot_data_scatter(MOFAobject, 
                  view = "IMC",
                  factor = 1,  
                  features = 4,
                  dot_size = 3,
                  sign = "positive",
                  color_by = "Group") + labs(y="Expression")

#This function generates a scatterplot of Factor 1 values (x-axis) versus expression values (y-axis) for the top 4 genes with largest negative weight. Samples are coloured by IGHV status:
plot_data_scatter(MOFAobject, 
                  view = "Luminex",
                  factor = 1,  
                  features = 4,
                  dot_size = 3,
                  sign = "negative",
                  color_by = "Group") + labs(y="Expression")

#The weights are useful to get an idea of which are top genes that drive the factors. 
#However, to get an idea of how well Factors are associated with genomic features 
#we can generate a scatter plot of the Factor values against mRNA expression for the genes with the highest weights:
#Top genes with positive weights
p <- plot_data_scatter(MOFAobject, 
                       view = "Luminex", 
                       factor = 1, 
                       features = 6,         # Number of features to show
                       sign = "positive",     # select top 6 features with positive weights
                       color_by = "Group",  # color cells by lineage
                       add_lm = T,          # add linear regression estimates
                       lm_per_group = F, 
                       dot_size = 2
)
p <- p + 
  scale_fill_manual(values=group.colors) +
  theme(legend.position = "none")
print(p)

#Top genes with negative weights
p <- plot_data_scatter(MOFAobject, 
                       view = "Luminex", 
                       factor = 1, 
                       features = 6,         # Number of features to show
                       sign = "negative",     # select top 6 features with negative weights
                       color_by = "Group",  # color cells by lineage
                       add_lm = T,          # add linear regression estimates
                       lm_per_group = F, 
                       dot_size = 2
)
p <- p + 
  scale_fill_manual(values=group.colors) +
  theme(legend.position = "none")
print(p)

#An alternative visualisation is to use a heatmap
plot_data_heatmap(MOFAobject, 
                  view = "Luminex",
                  factor = 1,  
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  annotation_samples = "Group",  
                  annotation_colors = list("Group"=group.colors), 
                  annotation_legend = F,
                  scale = "row")

#plot_data_heatmap has an interesting argument to “beautify” the heatmap: denoise = TRUE. 
#Instead of plotting the (noisy) input data, we can plot the data reconstructed by the model, where noise has been removed:
plot_data_heatmap(MOFAobject, 
                  view = "Luminex",
                  factor = 1,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  annotation_samples = "Group",  
                  annotation_colors = list("Group"=group.colors), 
                  annotation_legend = F,
                  scale = "row")

#Generate Heatmap of the top 25 features with largest loading, with the imputation of missing values(impute=TRUE):
df <- MOFAobject@samples_metadata[,c("sample","Group")]
rownames(df) <- df$sample; df$sample <- NULL

plot_data_heatmap(MOFAobject, 
                  view = "Luminex", 
                  factor = 2, 
                  impute = TRUE,
                  features = 25,
                  # extra arguments passed to `pheatmap`
                  denoise = TRUE,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  annotation_samples = "Group",  
                  annotation_colors = list("Group"=group.colors), 
                  annotation_legend = F,
                  scale = "row")

#4.4- Factor scatterplot
#A scatterplot of the first two factors shows the embedding of timepoint-species combination in the latent space. 
#Again we see, that the aligned times show a gradual development along the developmental trajectory captured by the conserved developmental programs on Factors 1 and 2.
plot_factors(MOFAobject, 2:3, color_by = "Group", shape_by = "Group", dot_size = 4) +
  scale_fill_manual(values = group.colors)

#We can also look at the factor values on the sample level. Here each dot correspond to one time-point combination.
gg_scatter <- plot_grid(
  plot_factors(MOFAobject, color_by = "Group") +
    theme(legend.position = "top") +
    scale_fill_manual(values = group.colors),
  plot_factors(MOFAobject, color_by = "Group") +
    theme(legend.position = "top") +
    scale_fill_manual(values = group.colors),
  plot_factors(MOFAobject, color_by = "Group") +
    theme(legend.position = "top"),
  nrow = 1, align = "h", axis = "tb")
gg_scatter

#4.4.1- Inspection of combinations of Factors
#Now that we have characterised the etiology of the two main Factors, let’s explore them together:
p <- plot_factors(MOFAobject, 
                  factors = c(2,3), 
                  color_by = "Group",
                  shape_by = "Group",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

#4.5- Interpolation
#Using the underlying Gaussian process for each factor we can interpolate to unseen time points for species that are missing data in these time points or intermediate time points. 
#If you want to show uncertainties, set only_mean to FALSE. 
#Note that these are only available if the new values for interpolation have been specified before training.
plot_interpolation_vs_covariate(MOFAobject_new, only_mean = FALSE)

plot_interpolation_vs_covariate(MOFAobject_new, only_mean = FALSE, factor = 1) +
  facet_wrap(~group)

# if the interpolation was not specified in mefisto_options before training
# or you want to interpolate to new values you can use interpolate_factors
# this however does not provide uncertainties
MOFAobject_new <- interpolate_factors(MOFAobject, new_values = seq(1, 14, 0.1))

## We recommend doing interpolation from python where additionally uncertainties are provided for the interpolation.
## Warning in interpolate_factors(MOFAobject, new_values = seq(1, 14, 0.1)): Object
## already contains interpolated factor values, overwriting it.
plot_interpolation_vs_covariate(MOFAobject_new)


#4.6- Non-linear dimensionality reduction from the MOFA factors
#The latent space inferred by MOFA can replace the PCA embedding as input to algorithms that learn non-linear manifolds such as t-SNE or UMAP. 
#This can be very useful to identify cellular populations and reconstruct complex pseudotime trajectories. 
#The advantage of  MOFA is that (1) we use information from all available omics, and (2) we can characterise the Factors and remove the technical ones.

#4.6.1- Run t-SNE
set.seed(42)
MOFAobject <- run_tsne(MOFAobject)

plot_dimred(MOFAobject,
            method = "TSNE",
            color_by = "Group"
) + scale_fill_manual(values=group.colors)

#4.6.2- Run UMAP
factors <- 1:get_dimensions(MOFAobject)[["K"]]
factors <- factors[!factors%in%c(4,7)]

mofa <- run_umap(MOFAobject, 
                 factors = factors, 
                 n_neighbors = 15,  
                 min_dist = 0.30)

plot_dimred(MOFAobject, 
            method = "UMAP", 
            color_by = "Group", 
            label = TRUE, 
            stroke=0.05, 
            dot_size = 5, 
            legend = FALSE
) + scale_fill_manual(values=group.colors)

#We can try to add some interpretatibility on the UMAP by visualising the contribution of each Factor on the different groups of cells.
for (i in paste0("Factor",1:3)) {
  p <- plot_dimred(MOFAobject, 
                   method = "UMAP", 
                   color_by = i, 
                   stroke = 0.05, 
                   dot_size = 1
  )
  print(p)
}


#4.7- Association with clinical covariates
#Factor 1 is correlated to pretty much all antibiotics, which is consistent with the observation that it separates the “Healthy non-antibiotic” individuals from all the others.
#Interestingly, Factor 4 is strongly correlated to Cephalosporin exposure.
correlate_factors_with_covariates(MOFAobject, 
                                  covariates = "DDOUD",
                                  plot = "r",  
)

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = "DDOUD",
                                  plot = "log_pval", # use "log_pval" to plot log p-values 
)

plot_factors(MOFAobject, 
             factors = c(1,2), 
             color_by = "TPO", 
             shape_by = "Group", 
             dot_size = 3.5
) + scale_fill_gradient(low = "white", high = "#BF3EFF")


#4.8- Building predictive models of clinical outcome (add the DDOUD to the model)
#4.8.1-Fit COX Models
library(survival)
library(survminer)

SurvObject <- Surv(MOFAobject@samples_metadata$DDOUD)
Z <- get_factors(MOFAobject)[[1]]
fit <- coxph(SurvObject ~ Z) 
fit

#4.8.1- Plot Hazard ratios
s <- summary(fit)
coef <- s[["coefficients"]]
p <- 

df <- data.frame(
  factor = factor(rownames(coef), levels = rev(rownames(coef))),
  p      = coef[,"Pr(>|z|)"], 
  coef   = coef[,"exp(coef)"], 
  lower  = s[["conf.int"]][,"lower .95"], 
  higher = s[["conf.int"]][,"upper .95"]
)

library(viridis)
library(dplyr)
mid<-0.05
highlight_df <- df %>% 
  filter(p<=0.01)

ggplot(df, aes(x=factor, y=coef, ymin=lower, ymax=higher)) +
  geom_pointrange(aes(color ="black")) + geom_point(aes(color ="black")) +
  geom_point(data=highlight_df, (aes (x=factor, y=coef, ymin=lower, ymax=higher, color ="red"))) +
  geom_pointrange(data=highlight_df, (aes (x=factor, y=coef, ymin=lower, ymax=higher, color ="red"))) +
  coord_flip() +
  scale_x_discrete() + 
  labs(y="Hazard Ratio - Death", x="") + 
  geom_hline(aes(yintercept=1), linetype="dotted") +
  theme_bw(base_size = 24) + scale_color_manual(values = c("black", "red"))

#scale_color_gradient2(midpoint=0.3, low="red", mid="white", high="blue", space ="RGB")
    
#4.8.2- Kaplan-Meier plots
#For illustration purposes we can also split the samples based on the factor values into two groups  
#using the maximally selected rank statistics from the maxstat R package and plot the Kaplan Meier plots per group.
library(maxstat)
df <- data.frame(
  time = SurvObject[,1], 
  event = SurvObject[,2], Z1 = Z[,2] #choose the column in Z related to the factor that has significant association with DDOUD
)
cut <- surv_cutpoint(df, variables='Z1')
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint

#The survfit function creates survival curves based on a formula:
fit2 <- survfit(Surv(time, event) ~ FactorCluster, df)
plot(survfit(Surv(time, event) ~ FactorCluster, df), 
     xlab = "Days of disease", 
     ylab = "Overall survival probability")

ggsurvplot(fit2, data = df,
           conf.int = TRUE, pval = TRUE,
           fun = function(y) y * 100, #linetype = "strata", risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "top", legend.labs = c(paste("Low IS 2"), paste("High IS 2")),
           xlab = "Days of disease", ylab="Survival probability (%)", title= "Factor 2")$plot + scale_color_manual(values = c("light blue", "salmon")) +
  scale_fill_manual(values = c("light blue", "salmon")) + scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  theme_classic(base_size = 22) + theme(legend.justification=c(),
                                   legend.position='bottom')

#Use summary with the times argument to check the survival probability levels
summary(survfit(Surv(time, event) ~ FactorCluster, df), times = 20)
#We get the log-rank p-value using the survdiff function. 
#For example, we can test whether there was a difference in survival time according to factor 2 levels
survdiff(Surv(time, event) ~ FactorCluster, df)

coxph(Surv(time, event) ~ FactorCluster, df)
#We can see a tidy version of the output using the tidy function from the broom package:
library(broom)
broom::tidy(
    coxph(Surv(time, event) ~ FactorCluster, df), 
    exp = TRUE
  ) %>% 
  table()
#Or use tbl_regression from the gtsummary package
library(gtsummary)
coxph(Surv(time, event) ~ FactorCluster, df) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

#Transform survival curves: plot cumulative events and hazard function
# Plot cumulative events
ggsurvplot(fit2, data = df,
           conf.int = TRUE, pval = FALSE,
           fun = "event", #linetype = "strata", risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "top", legend.labs = c(paste("Low Factor 2"), paste("High Factor 2")),
           xlab = "Days of disease", ylab="Cumulative events (death)", title= "Factor 2")$plot + scale_color_manual(values = c("light blue", "salmon")) +
  scale_fill_manual(values = c("light blue", "salmon")) + scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12)

#Plot the cumulative hazard function
ggsurvplot(fit2, data = df,
           conf.int = TRUE, pval = TRUE,
           fun = "cumhaz", #linetype = "strata", risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "top", legend.labs = c(paste("Low Factor 2"), paste("High Factor 2")),
           xlab = "Days of disease", ylab="cumulative hazard", title= "Factor 2")$plot + scale_color_manual(values = c("light blue", "salmon")) +
  scale_fill_manual(values = c("light blue", "salmon")) + scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12)


#Survival analysis with Factor 1
df <- data.frame(
  time = SurvObject[,1], 
  event = SurvObject[,2], Z1 = Z[,1] #choose the column in Z related to the factor that has significant association with DDOUD
)
cut <- surv_cutpoint(df, variables='Z1')
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint

#The survfit function creates survival curves based on a formula:
fit3 <- survfit(Surv(time, event) ~ FactorCluster, df)
plot(survfit(Surv(time, event) ~ FactorCluster, df), 
     xlab = "Days of disease", 
     ylab = "Overall survival probability")

ggsurvplot(fit3, data = df,
           conf.int = TRUE, pval = TRUE,
           fun = function(y) y * 100, #linetype = "strata", risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "top", legend.labs = c(paste("High Factor 1"), paste("Low Factor 1")),
           xlab = "Days of disease", ylab="Survival probability (%)", title= "Factor 1")$plot + scale_color_manual(values = c("light blue", "salmon")) +
  scale_fill_manual(values = c("light blue", "salmon")) + scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12)

#Survival analysis with Factor 3
df <- data.frame(
  time = SurvObject[,1], 
  event = SurvObject[,2], Z1 = Z[,3] #choose the column in Z related to the factor that has significant association with DDOUD
)
cut <- surv_cutpoint(df, variables='Z1')
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint

#The survfit function creates survival curves based on a formula:
fit4 <- survfit(Surv(time, event) ~ FactorCluster, df)
plot(survfit(Surv(time, event) ~ FactorCluster, df), 
     xlab = "Days of disease", 
     ylab = "Overall survival probability")

ggsurvplot(fit4, data = df,
           conf.int = TRUE, pval = TRUE,
           fun = function(y) y * 100, #linetype = "strata", risk.table = TRUE, risk.table.y.text.col = TRUE,
           legend = "top", legend.labs = c(paste("Low Factor 3"), paste("High Factor 3")),
           xlab = "Days of disease", ylab="Survival probability (%)", title= "Factor 3")$plot + scale_color_manual(values = c("light blue", "salmon")) +
  scale_fill_manual(values = c("light blue", "salmon")) + scale_x_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12)

#5- Smoothness and sharedness of factors: Investigating the learnt hyper parameters
#In addition to the factor values and the alignment the model also inferred an underlying Gaussian process that generated these values.
#By looking into it we can extract information on the smoothness of each factor, i.e. how smoothly it varies along developmental time, as well as the sharedness of each factor, 
#i.e. how much the species (groups) show the same underlying developmental pattern and how the shape of their developmental trajectory related to a given developmental module (Factor) clusters between species.

#5.1 Smoothness score
#The scale parameters of the Gaussian process capture the smoothness of the model. A scale of 1 indicates high smoothness, a scale of 0 variation independent of time. Here all factors show a very high degree of smoothness.
get_scales(MOFAobject)

#For visualization we visualize the scale of the Gaussian process as a smoothness score in a barplot for each factor. 
#No color indicates no smoothness, a fully colored bar a very smooth factor.
plot_smoothness(MOFAobject) 

#5.2- Sharedness
#The group kernel of the Gaussian process can give us insights into the extent to which a temporal pattern is shared across species for the developmental module captured by each process. 
#Here, we see a strong clustering of the first processes, indicating that these are conserved in evolution of the species, while the last two processes show distinct clusters.

# We can use the following, below we make some more custom plot
plot_group_kernel(MOFAobject, factors = "all",
                  cellheight =  15, cellwidth = 15,
                  treeheight_col = 3, treeheight_row = 3)

#For visualization we calculate an overall sharedness score per factor between 0 and 1. No color indicates no sharedness, a fully colored bar a fully shared factor.
plot_sharedness(MOFAobject)

#6- Alignment
#Next, we further look into the alignment that was learnt by the model by plotting the mapping of mouse developmental stages to the stages of other species.
plot_alignment(MOFAobject)

#7- Extracting data for downstream analysis
#The user can extract the feature weights, the data and the factors to generate their own plots.
#Extract factors
# "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
factors <- get_factors(MOFAobject, factors = "all")
lapply(factors,dim)
factors3 <- factors[["single_group"]]
write.csv(factors3,'MOFA_MEFISTO_factor_values.csv')

#Extract weights
# "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
weights <- get_weights(model, views = "all", factors = "all")
lapply(weights,dim)

#Extract data
# "data" is a nested list of matrices, one matrix per view and group with dimensions (nfeatures, nsamples)
data <- get_data(model)
lapply(data, function(x) lapply(x, dim))[[1]]

#For convenience, the user can extract the data in long data.frame format:
factors3 <- get_factors(MOFAobject, as.data.frame = T)
head(factors, n=3)

weights <- get_weights(model, as.data.frame = T)
head(weights, n=3)

data <- get_data(model, as.data.frame = T)
head(data, n=3)



#8- Means and Standard deviation calculation for each factor 
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
multiplex <- read.csv("MOFA_MEFISTO_factor_values.csv")

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


df1 <- data_summary(multiplex, varname="Factor1", 
                    groupnames=c("Group", "Day_of_Sampling"))
head(df1)

df2 <- data_summary(multiplex, varname="Factor2", 
                    groupnames=c("Group", "Day_of_Sampling"))
head(df2)

df3 <- data_summary(multiplex, varname="Factor3", 
                    groupnames=c("Group", "Day_of_Sampling"))
head(df3)


#Standard error or mean calculation for each factor 
df1_1 <- summarySE(multiplex, measurevar="Factor1", groupvars=c("Group", "Day_of_Sampling"))
head(df1_1)

df2_1 <- summarySE(multiplex, measurevar="Factor2", groupvars=c("Group", "Day_of_Sampling"))
head(df2_1)

df3_1 <- summarySE(multiplex, measurevar="Factor3", groupvars=c("Group", "Day_of_Sampling"))
head(df3_1)


write.csv(df1_1, "MOFA_MEFISTO_Factor1_means.csv")
write.csv(df2_1, "MOFA_MEFISTO_Factor2_means.csv")
write.csv(df3_1, "MOFA_MEFISTO_Factor3_means.csv")

#Statistical test
compare_means(Factor1 ~ Group, data = multiplex, method = "wilcox.test", p.adjust.method	= 'fdr', group.by = "Day_of_Sampling")

compare_means(Factor2 ~ Group, data = multiplex, method = "wilcox.test", p.adjust.method	= 'fdr', group.by = "Day_of_Sampling")

compare_means(Factor3 ~ Group, data = multiplex,  method = "wilcox.test", p.adjust.method	= 'fdr', group.by = "Day_of_Sampling")


