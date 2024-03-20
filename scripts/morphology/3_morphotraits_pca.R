# Analysis of morphological traits in _Aegilops_-_Triticum_ species
# Part 3. PCA
# author: Concetta Burgarella (Uppsala University)
# date: September 2021, updated November 2022

#### The objective of this script is to characterize the mating system of 13 _Aegilops_/_Triticum_ species using morpho-functional traits related to reproduction. We will use a set of traits that are expected to be more impacted by the rate of self-reproduction in anemophilous species.

#packages to install
if (!require("psych"))   install.packages("psych", dependencies = TRUE)
if (!require("pals"))   install.packages("pals", dependencies = TRUE)
if (!require("ade4"))   install.packages("ade4", dependencies = TRUE)
if (!require("dplyr"))   install.packages("dplyr", dependencies = TRUE)
if (!require("adegenet"))   install.packages("adegenet", dependencies = TRUE)
if (!require("MASS"))   install.packages("MASS", dependencies = TRUE)
if (!require("pals"))   install.packages("pals", dependencies = TRUE)
if (!require("ape"))   install.packages("ape", dependencies = TRUE)


# 1) Load and subset datasets to retain target species and variables ####
load("outputs/morphology/imputed_morpho_data.Rdata")


### Select species

#### We extract the data for the 13 _Aegilops_/_Triticum_ species for the analysis. \

# List of species included in the table
unique(imputed_data$species_code)
# Create the data set *Aegilops/Triticum* (the first 93 individuals)
aeg <- imputed_data[c(grep("A", imputed_data$species_code), grep("Tmo", imputed_data$species_code), grep("Tur", imputed_data$species_code)),]
aeg$species_code <- as.character(aeg$species_code) # transform type of variable to "character"
table(aeg$species_code)
species <- unique(aeg$species_code) # list of species name (in alphabetical order)

### Select morphological variables 
# We will use 6 variables describing traits that we suppose to be most affected by species mating system: \

# Eta_len_norm: mean anther length normalised by the flower length. \
# Eta_wid_norm: mean anther width normalised by the flower length. \
# Sti_norm: mean stigma length normalised by the flower length. \
# maleinvest: male investiment, calculated as mean anther length/mean ovary length. \
# Compactness: a measure of spikelet overlap calculated as the ratio of (mean flower length* number_fertile_flower)/mean spikelet length. \
# AF_index: autonomous seed set, calculated as self-fertilised_seed_number/(self-fertilized_spikelet_number*number_fertile_flower/spikelet). \

# List of variables included in the table
colnames(imputed_data)[-(1:3)]
# Subset variables 
attach(aeg)
var <- as.data.frame(cbind(eta_len_norm,eta_wid_norm,sti_norm,compactness,AF_index,maleinvest))
detach(aeg)
# Look at pairwise correlation between traits 
library(PerformanceAnalytics)
chart.Correlation(var, histogram = TRUE, method = "pearson")

#### Set a palette of colors (palette1)
# We want to use a palette that aims to have at its extremes bluish colors (meant for outcrossing species *A. speltoides* and *A. mutica*) and reddish colors (meant for selfing ones *A. tauschii*, *T. urartu*). In the middle yellow/green colors for mixed ones.
pal.bands(kovesi.diverging_rainbow_bgymr_45_85_c67(13))


# 2) PCA of reproductive variables ####

pca_individuals <- dudi.pca(var, scale=T, scannf = F, nf = 4)                # n axes = 2
inertie <- pca_individuals$eig/sum(pca_individuals$eig)*100

#### Assign colors to each species based on its position on PC1
# Create vector of colors for the 13 Aegilops/Triticum species
mycol <- rep(NA, nrow(aeg))
mycol[c(which(aeg$species_code == "Amu"))] <- "#085CF8"
mycol[c(which(aeg$species_code == "Asp"))] <- "#1684A6"   # "royalblue"
mycol[c(which(aeg$species_code == "Ash"))] <- "#3C9E49"
mycol[c(which(aeg$species_code == "Aca"))] <- "#65AF1E"
mycol[c(which(aeg$species_code == "Alo"))] <- "#98BB18"
mycol[c(which(aeg$species_code == "Aum"))] <- "#C7C612" 
mycol[c(which(aeg$species_code == "Abi"))] <- "#F3CC1D"
mycol[c(which(aeg$species_code == "Aco"))] <- "#FDAF55"
mycol[c(which(aeg$species_code == "Tmo"))] <- "#FE8F7B"
mycol[c(which(aeg$species_code == "Ase"))] <- "#FC6A9B"
mycol[c(which(aeg$species_code == "Aun"))] <- "#F64497" 
mycol[c(which(aeg$species_code == "Ata"))] <- "#E92952"
mycol[c(which(aeg$species_code == "Tur"))] <- "#D70500"

# Plot PCA
op <- par(mfrow=c(1,2), cex=0.8)
# 1. plot space of variables
s.corcircle(pca_individuals$co)  
# 2. scatterplot PC1-PC2
x= pca_individuals$li[,1]
y= pca_individuals$li[,2]
plot(x, y, col=mycol, cex=0.5, pch=19, xlab=paste("PC1 (",round(inertie[1]),"%)",sep=""), ylab=paste("PC2 (",round(inertie[2]),"%)",sep=""))
# Add centroid position and species name
centroid_x <- tapply(x, aeg$species_code,mean)
centroid_y <- tapply(y, aeg$species_code,mean)
points(centroid_x, centroid_y, pch=18, col=mycol[match(species,aeg$species_code)], cex=1.6)
# tune the y position
centroid_y_text <- centroid_y + 0.25
centroid_x_text <- centroid_x - 0.2
centroid_x_text[11:13] <- centroid_x[11:13] + 0.2
text(centroid_x_text,centroid_y_text,labels=species, col=mycol[match(species,aeg$species_code)], cex=1, font=2)
par(op)

pdf(file = paste("figures/main/PCA_pc1_pc2_final_palette.pdf",sep=""), pointsize = 15, height=7, width=10)
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(0.85, 0.15))
# plot 2. scatterplot PC1-PC2
op <- par(mai=c(1,1.2,0.3,0.006))
x= pca_individuals$li[,1]
y= pca_individuals$li[,2]
plot(x, y, col=mycol, cex.lab=1.4, cex=0.6, pch=19, #pch=patterns,
     xlab=paste("PC1 (",round(inertie[1]),"%)",sep=""), ylab=paste("PC2 (",round(inertie[2]),"%)",sep=""))
# Add centroid position
centroid_x <- tapply(x, aeg$species,mean)
centroid_y <- tapply(y, aeg$species,mean)
points(centroid_x, centroid_y, pch=18, col=mycol[match(species,aeg$species)], cex=2.6)
par(op)
# plot 3. legend
par(mai=c(0.08,0.01,0.4,0.01))
plot(x, y, type="n",bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
legend(-4.4,4, bty="n", legend=species, col = mycol[match(species,aeg$species)], pch=18, pt.cex=2,cex=1.4)
par(op)
dev.off()

pdf(file = paste("figures/sup_mat/Suppl_Figure_Sx_PCA_variables_space.pdf"), pointsize=17, height=9, width=10)
s.corcircle(pca_individuals$co)
dev.off()


# 3) Correlation between principal component PC1 and fixation index _F_ ####

#### We now want to quantify the correlation between the selfing syndrome described by PC1 and species fixation index _F_ with a linear regression analysis.
# Load polymorphism values
polym <- as.data.frame(t(read.delim("data/polymorphism/Summary_polymorphism_Triticeae_extended.csv",sep = ";", row.names=1)))
# subset the table to retain only Aegilops/triticum species 
polym_aeg <- polym[match(species,row.names(polym)),]
# Perform the linear regression 
pc1 <- pca_individuals$li[,1]
pc1_mean <- tapply(pc1,aeg$species_code, mean)
x <- pc1_mean
y <- polym_aeg$`Mean_Fit_W&C`
regression <- lm(y ~ x)
summary(regression)

# Plot regression F on PC1
plot(x,y, xlab="PC1",ylab=substitute(paste(italic("F"))), pch=18, xlim=c(min(x)-0.15, max(x)+0.55), ylim = c(0.18, 1.05), col=mycol[match(species,aeg$species_code)], cex=1.4)
abline(regression, lty=2 )
R2 <- summary(regression)$r.squared
text(x=x+0.3, y=y, labels=species, cex=1, col=mycol[match(species,aeg$species_code)] )
legend(x=-1,y=0.4, legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.2)
#### The correlation is very strong (R-squared:  0.8262, p-value: 1.686e-05).
#### Since the different species may have different distributions of _F_, we repeat teh regression by applying a weight to each _F_ value that is inversely proportional to the variance of the estimation. \ 

### Weighting regression by the variance of _F_ \
# Assign a weight to each species estimate of _F_ that depends on its variance.
# Calculate variance from CI for F
N <- polym_aeg$Sample_size_for_Fit
CIupper <- polym_aeg$Fit_CI_upper
CIlower <- polym_aeg$Fit_CI_lower
variance <- sqrt(N)*(CIupper-CIlower)/3.92
# Perform the linear regression 
pc1 <- pca_individuals$li[,1]
pc1_mean <- tapply(pc1,aeg$species_code, mean)
x <- pc1_mean
y <- polym_aeg$`Mean_Fit_W&C`
regression <- lm(y ~ x, weights = 1/variance)
summary(regression)
#### After assigning a weight to _F_, the correlation is still very strong and significant, R-squared:  0.7727, p-value: 7.566e-05. \

#### Another important factor that could drive the strong correlation is the phylogenetic relationships among species. Thus, we check if the correlation is still significant after correcting for it. \ 
### Applying phylogenetically independent contrast correction \
# Ultrametric tree retrieved from Glémin et al. 2019
# Note that the individual names have been homogeneized: a different code was used for boeoticum with TSxx instead of Trxx.
# Int he input tree file TS was replaced by Tr
TREEFILE <- "data/MLtree_OneCopyGenes_ultrametric.tree"
fulltree <- read.tree(TREEFILE)
# Extraction of one individual per species by hand
tokeep <- fulltree$tip.label[c(5,9,13,17,21,25,28,30,33,36,40,42,45)]
tree <- keep.tip(fulltree,tokeep)
# Remove individual code in species name
newlabels <- sapply(tree$tip.label,function(x) strsplit(x,"_Tr")[[1]][1])
tree$tip.label <- newlabels
plot(tree)

# Function to apply phylogenetically independent contrast directly
# The argument of this function is a numeric vector of size 13 with species names as label for each value
# It returns a 12 row, two-column matrix, with the phylogenetically independent contrasts in the first column and their expected variance in the second column
PIC <- function(x) pic(x,tree, var.contrasts = T)
# Application
# The PIC function is applied to the two vectors to be compared, then analyses can be done with the contrast values instead of the raw values.
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x2 <- PIC(x)
y2 <- PIC(y)
plot(x2,y2, col="blue")
abline(lm(y2~x2))
# Run the linear regression model with the correction for phylogeny
summary(lm(y2[,1] ~ x2[,1]))
#### Correcting for phylogenetic relatedness the correlation remains strong and significant (R-squared:  0.7753, p-value: 0.0001565). \

#### It is also possible to correct the regression with contrasts taking into account the variance of the contrast estimate. \
### Phylogenetically independent contrast + correction for contrast variance
# Run the linear regression model with the correction for phylogeny weighting the contrast estimate of y by its variance
summary(lm(y2[,1] ~ x2[,1], weights = 1/y2[,2]))
#### The correlation is still very strong (R2=0.707) and significant (p-value: 0.0006118).


#### Plot with confidence intervals for _F_ : not nice figure, not published
library("plotrix")
pc1 <- pca_individuals$li[,1]
pc1_mean <- tapply(pc1,aeg$species_code, mean)
x <- pc1_mean
y <- polym_aeg$`Mean_Fit_W&C`
plot(x,y, xlab="PC1",type="n", ylab=substitute(paste(italic("F"))), xlim=c(min(x)-0.15, max(x)+0.55), ylim = c(0.18, 1.05))
plotCI(x,y, li=polym_aeg$Fit_CI_lower, ui=polym_aeg$Fit_CI_upper, add=T, col=mycol[match(species,aeg$species_code)], cex=1)
abline(lm(y ~ x), lty=2 )
R2 <- summary(lm(y ~ x))$r.squared
text(x=x+0.3, y=y, labels=species, cex=1, col=mycol[match(species,aeg$species_code)] )
legend(x=0,y=0.4, legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.2)

# Figure export
pdf(paste0("figures/main/Fis_pc1_final_palette.pdf",sep=""), pointsize = 20)
pc1 <- pca_individuals$li[,1]
pc1_mean <- tapply(pc1,aeg$species_code, mean)
x <- pc1_mean
y <- polym_aeg$`Mean_Fit_W&C`
plot(x,y, xlab="PC1",ylab=substitute(paste(italic("F"))), pch=18, xlim=c(min(x)-0.15, max(x)+0.15), ylim = c(0.18, 1.05), col=na.omit(mycol[match(species,aeg$species)]), cex=2.2)
abline(lm(y ~ x), lty=2 )
R2 <- summary(lm(y ~ x))$r.squared
legend(x=-1.5,y=0.35, legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n", cex=1)
legend(x=-2.5,y=0.25, legend = bquote("(With phylogenetic correction"~R^2 == "0.78***)"), bty = "n", cex=0.6) # with phylogenetically independent contrast correction
dev.off()

# Save the PCA results for use in subsequent analyses 
save(list=c("pca_individuals", "aeg"), file="outputs/morphology/PCA_morphological_traits.Rdata")
