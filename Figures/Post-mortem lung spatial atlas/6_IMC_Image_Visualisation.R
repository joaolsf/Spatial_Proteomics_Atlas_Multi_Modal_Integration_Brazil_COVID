library(imcRtools)
library(cytomapper)
library(openxlsx)
library(stringr)
library(dittoSeq)
library(RColorBrewer)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)
library(ggplot2)
library(scran)
library(CATALYST)
library(kohonen)
library(ConsensusClusterPlus)
library(patchwork)
library(pheatmap)
library(gridExtra)
library(SingleCellExperiment)
library(tidyverse)
library(ggridges)

#1- Load data
#First, we will read in the previously generated SpatialExperiment object.

#The following section describes how to visualize the abundance of biomolecules (e.g. protein or RNA) 
#as well as cell-specific metadata on images. 
#Section 11.1 focuses on visualizing pixel-level information including the generation of pseudo-color composite images. 
#Section 11.2 highlights the visualization of cell metadata (e.g. cell phenotype) as well as summarized pixel intensities on cell segmentation masks.

#The cytomapper R/Bioconductor package was developed to support the handling and visualization of multiple multi-channel images and segmentation masks (Eling et al. 2020). 
#The main data object for image handling is the CytoImageList container which we used in Section 5 to store multi-channel images and segmentation masks.
#We will first read in the previously processed data and randomly select 3 images for visualization purposes.

library(SpatialExperiment)
library(cytomapper)
cytomapper_sce <- readRDS("cytomapper_sce.rds")
images <- readRDS("images.rds")
masks <- readRDS("masks.rds")

# Sample images
set.seed(220517)
cur_id <- sample(unique(cytomapper_sce$sample_id), 3)
cur_images <- images[names(images) %in% cur_id]
cur_masks <- masks[names(masks) %in% cur_id]

#1.1- Pixel visualization

#The following section gives examples for visualizing individual channels or multiple channels as pseudo-color composite images. 
#For this the cytomapper package exports the plotPixels function which expects a CytoImageList object storing one or multiple multi-channel images. 
#In the simplest use case, a single channel can be visualized as follows:

#The bcg parameter (default c(0, 1, 1)) stands for “background”, “contrast”, “gamma” and controls these attributes of the image. 
#This parameter takes a named list where each entry specifies these attributes per channel. 
#The first value of the numeric vector will be added to the pixel intensities (background); 
#pixel intensities will be multiplied by the second entry of the vector (contrast); 
#pixel intensities will be exponentiated by the third entry of the vector (gamma). 
#In most cases, it is sufficient to adjust the second (contrast) entry of the vector.
cur_images2 <- images[names(images) %in% "C27_ROI1"]
cur_masks2 <- masks[names(masks) %in% "C27_ROI1"]
                   
plotPixels(cur_images, 
           colour_by = "SARSCoV2",
           bcg = list(SARSCoV2 = c(0, 5, 1)))

#The following example highlights the visualization of 6 markers (maximum allowed number of markers) at once per image. 
#The markers indicate the spatial distribution of tumor cells (E-caherin), T cells (CD3), B cells (CD20), CD8+ T cells (CD8a),
#plasma cells (CD38) and proliferating cells (Ki67).

plotPixels(cur_images, 
           colour_by = c("SARSCoV2", "PanCK", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      PanCk = c(0, 5, 1),
                      DNA = c(0, 5, 1)))
           
#CD8a = c(0, 5, 1),
#CD38 = c(0, 8, 1),
#Ki67 = c(0, 5, 1))

#1.1.1- Adjusting colors
#The default colors for visualization are chosen by the additive RGB (red, green, blue) color model. 
#For six markers the default colors are: red, green, blue, cyan (green + blue), magenta (red + blue), yellow (green + red). 
#These colors are the easiest to distinguish by eye. 
#However, you can select other colors for each channel by setting the colour parameter:

#The colour parameter takes a named list in which each entry specifies the colors from which a color gradient is constructed via colorRampPalette. 
#These are usually vectors of length 2 in which the first entry is "black" and the second entry specifies the color of choice. 
#Although not recommended, you can also specify more than two colors to generate a more complex color gradient.

plotPixels(cur_images, 
           colour_by = c("SARSCoV2", "PanCK", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      PanCK = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "burlywood1"),
                         PanCK = c("black", "cyan2"),
                         DNA1 = c("black", "firebrick1")))

plotPixels(cur_images, 
           colour_by = c("SARSCoV2", "CD206", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      CD206 = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "yellow"),
                         CD206 = c("black", "firebrick1"),
                         DNA1 = c("black", "cyan2")))


#1.1.2- Image normalization
#As an alternative to setting the bcg parameter, images can first be normalized. 
#Normalization here means to scale the pixel intensities per channel between 0 and 1 (or a range specified by the ft parameter in the normalize function). 
#By default, the normalize function scales pixel intensities across all images contained in the CytoImageList object (separateImages = FALSE). 
#Each individual channel is scaled independently (separateChannels = TRUE).

#After 0-1 normalization, maximum pixel intensities can be clipped to enhance the contrast of the image (setting the inputRange parameter). 
#In the following example, the clipping to 0 and 0.2 is the same as multiplying the pixel intensities by a factor of 5.

# 0 - 1 channel scaling across all images
ED_image <- images[names(images) %in% c("C27_ROI1", "C27_ROI2", "C27_ROI3", "C27_ROI4")]
ED_masks <- masks[names(masks) %in% c("C27_ROI1", "C27_ROI2", "C27_ROI3", "C27_ROI4")]

LD_image <- images[names(images) %in% c("C16_ROI1", "C16_ROI2", "C16_ROI3", "C16_ROI4")]
LD_masks <- masks[names(masks) %in% c("C16_ROI1", "C16_ROI2", "C16_ROI3", "C16_ROI4")]

plotPixels(ED_image, 
           colour_by = c("SARSCoV2", "CD206", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      CD206 = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "yellow"),
                         CD206 = c("black", "firebrick1"),
                         DNA1 = c("black", "cyan2")))

plotPixels(ED_image, 
           colour_by = c("SARSCoV2", "CD163", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      CD163 = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "yellow"),
                         CD163 = c("black", "firebrick1"),
                         DNA1 = c("black", "cyan2")))

plotPixels(LD_image, 
           colour_by = c("SARSCoV2", "PanCK", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      PanCK = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "yellow"),
                         PanCK = c("black", "firebrick1"),
                         DNA1 = c("black", "cyan2")))


# Clip channel at 0.2
norm_ED_image <- normalize(ED_image)
norm_ED_image2 <- normalize(norm_ED_image, inputRange = c(0, 0.2))

plotPixels(norm_ED_image, 
           colour_by = c("SARSCoV2", "CD163", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      CD163 = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "yellow"),
                         CD163 = c("black", "firebrick1"),
                         DNA1 = c("black", "cyan2")))

plotPixels(norm_ED_image, 
           colour_by = c("SARSCoV2", "CD206", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      CD206 = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "yellow"),
                         CD206 = c("black", "firebrick1"),
                         DNA1 = c("black", "cyan2")))

norm_LD_image <- normalize(LD_image)
norm_LD_image2 <- normalize(norm_ED_image, inputRange = c(0, 0.2))

plotPixels(norm_LD_image, 
           colour_by = c("SARSCoV2", "PanCK", "DNA1"),
           bcg = list(SARSCoV2 = c(0, 5, 1),
                      PanCK = c(0, 5, 1),
                      DNA1 = c(0, 5, 1)),
           colour = list(SARSCoV2 = c("black", "yellow"),
                         PanCK = c("black", "firebrick1"),
                         DNA1 = c("black", "cyan2")))










