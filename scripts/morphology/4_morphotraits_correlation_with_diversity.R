# Analysis of morphological traits in _Aegilops_-_Triticum_ species
# Part 4. Correlation between morphological traits and polymorphism patterns
# Concetta Burgarella (Uppsala University, Sweden) and Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022

#### The objective of this script is to analyze the correlation between the phenotypic selfing syndrome (described by a PCA of reproductive traits) and genetic diversity, and compare with species range as driver of diversity.

#packages to install
if (!require("psych"))   install.packages("psych", dependencies = TRUE)
if (!require("pals"))   install.packages("pals", dependencies = TRUE)
if (!require("ade4"))   install.packages("ade4", dependencies = TRUE)
if (!require("dplyr"))   install.packages("dplyr", dependencies = TRUE)
if (!require("adegenet"))   install.packages("adegenet", dependencies = TRUE)
if (!require("MASS"))   install.packages("MASS", dependencies = TRUE)
if (!require("plotrix"))   install.packages("plotrix", dependencies = TRUE)
if (!require("ape"))   install.packages("ape", dependencies = TRUE)


# 1) Prepare datasets ####

# Load polymorphism values
polym <- as.data.frame(t(read.delim("data/polymorphism/Summary_polymorphism_Triticeae_extended.csv", row.names=1, sep=";")))
# Load Principal Component Analysis of reproductive traits
load("outputs/morphology/PCA_morphological_traits.Rdata")

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

# Make a list of species name (in alphabetical order)
species <- unique(aeg$species_code) 


# 2) Correlation with polymorphism statistics ####

#### We use the first principal component (PC1) summarising species diversity in reproductive traits as a morphological descriptor of the selfing syndrome and check how polymorphism levels correlate with it.
### **Synonymous polymorphism _$\pi_S$_**
#### We quantify first how the reproductive traits explain the species wide neutral nucleotide diversity, by looking at the correlation between synonymous polymorphism _$\pi_S$_ and PC1
pc1 <- pca_individuals$li[,1]
pc1_mean <- tapply(pc1,aeg$species_code, mean)
corr <-  match(species,row.names(polym))

layout(matrix(c(1,2),ncol=2))
# piS
x <- pc1_mean
y <- polym$Mean_piS[corr]
plot(x,y, xlab="PC1",ylab=expression(italic(pi)[S]), pch=18, xlim=c(min(x)-0.15, max(x)+0.15), 
     ylim = c(min(y), max(y)+0.05*max(y)), col=mycol[match(species,aeg$species_code)], cex=1.4)
# Add a weight to the regression (inverse of the variance of the estimator of piS)
# w <- 1/ (polym$piS_CI_upper[corr] -polym$piS_CI_lower[corr])^2
w <- NULL # unweighted regression
abline(lm(y ~ x,weights = w), lty=2 )
R2 <- summary(lm(y ~ x,weights = w))$r.squared
text(x=x, y=y+0.0008, labels=species, cex=0.8,col=mycol[match(species,aeg$species_code)])
legend(max(x)-2,y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.2)


# Run the linear regression model (with piS)
summary(lm(y ~ x,weights = w))

# log10(piS)
x <- pc1_mean
y <- log10(polym$Mean_piS[corr])
plot(x,y, type="n",  xlab="PC1",ylab=expression(italic(pi)[S]~(log[10]~scale)), xlim= c(min(x)-0.15, max(x)+0.15), 
     ylim = c(min(y),max(y)-0.05*max(y)), pch=18, cex=1.4, yaxt="n", col=mycol[match(species,aeg$species_code)], ) #, yaxt="n")
plotCI(x,y, li=log10(polym$piS_CI_lower[corr]), ui=log10(polym$piS_CI_upper[corr]), add=T,
       col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, lwd=2)
axis(side=2, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05), cex.axis=0.8,padj=0)

# Add a weight to the regression (inverse of the variance of the estimator of log(piS): var(piS)/mean(piS)^2
# w <- polym$Mean_piS[corr]/ (polym$piS_CI_upper[corr] -polym$piS_CI_lower[corr])^2
w <- NULL # unweighted regression
abline(lm(y ~ x,weights = w), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x,weights = w))$r.squared
text(x=x[-c(2:5,11)] - 0.35, y=y[-c(2:5,11)], labels=species[-c(2:5,11)], cex=1,
     col=mycol[match(species[-c(2:5,11)],aeg$species_code)]) # labels on the left
text(x=x[c(2:5,11)], y=y[c(2:5,11)] + 0.07, labels=species[c(2:5,11)], cex=1,
     col=mycol[match(species[c(2:5,11)],aeg$species_code)]) # labels on top
legend(max(x)-2,y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.2)

# Run the linear regression model (with piS)
summary(lm(y ~ x,weights = w))


#### Now we run the same models, but correcting for the phylogenetic relationships among species.
# Ultrametric tree retrieved from Glémin et al. 2019
# Note that the individual names have been homogenized: a different code was used for boeoticum with TSxx instead of Trxx.
# Int he input tree file TS was replaced by Tr
fulltree <- read.tree("data/MLtree_OneCopyGenes_ultrametric.tree")
#plot(fulltree)

# Extraction of one individual per species
# List of tips
fulltree$tip.label
# Extraction by hand
tokeep <- fulltree$tip.label[c(5,9,13,17,21,25,28,30,33,36,40,42,45)]
tree <- keep.tip(fulltree,tokeep)
# Remove individual code in species name
newlabels <- sapply(tree$tip.label,function(x) strsplit(x,"_Tr")[[1]][1])
tree$tip.label <- newlabels
# plot(tree)

# The pic function is applied to the two vectors to be compared, then analyses can be done with the contrast values instead of the raw values.
# The analyse is done on the log(piS)
x <- pc1_mean
y <- log10(polym$Mean_piS[corr])
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(y,tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
# w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
plot(x.pic,y.pic)
abline(lm(y.pic~x.pic-1,weights = w.pic))
summary(lm(y.pic~x.pic-1, weights = w.pic))
#### Synonymous nucleotide diversity is strongly correlated with the morphological reproductive features Adjusted R-squared: 0.75, p-value:0.00014), with species showing the strongest selfing syndrome being characterized by lower genetic diversity at the genome scale. Significant correlation stands even after correcting for the phylogenetic relationships among species (Adjusted R-squared:0.67, p-value: 0.00064).
#### This relationship doesn't change if the median of PC1 per species is used instead of the mean (only the two most extreme species in the PC1 axis slightly change position - *Ae. mutica* and *Ae. tauschii* - but the global picture remains the same). Not shown

pc1 <- pca_individuals$li[,1]
pc1_median <- tapply(pc1,aeg$species_code, median)
corr <-  match(species,row.names(polym))
layout(matrix(c(1,2),ncol=2))
# piS
x <- pc1_median
y <- polym$Mean_piS[corr]
plot(x,y, xlab="PC1",ylab=expression(italic(pi)[S]), pch=18, cex=1.4, xlim=c(min(x)-0.15, max(x)+0.15), ylim = c(min(y), max(y)+0.05*max(y)), col=mycol[match(species,aeg$species_code)])
abline(lm(y ~ x), lty=2 )
R2 <- summary(lm(y ~ x))$r.squared
text(x=x, y=y+0.0008, labels=species, cex=0.8,col=mycol[match(species,aeg$species_code)])
legend(max(x)-2,y=max(y), legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n")
summary(lm(y ~ x))

# log10(piS)
x <- pc1_median
y <- log10(polym$Mean_piS[corr])
plot(x,y, xlab="PC1",ylab=expression(italic(pi)[S]~(log[10]~scale)), xlim= c(min(x)-0.15, max(x)+0.15), ylim = c(min(y), max(y)-0.05*max(y)), pch=18, cex=1.4, yaxt="n", col=mycol[match(species,aeg$species_code)])
axis(side=2, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05), las=2, cex.axis=0.8)
abline(lm(y ~ x), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x))$r.squared
text(x=x, y=y+0.05, labels=species, cex=0.8,col=mycol[match(species,aeg$species_code)]) 
legend(max(x)-2,y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n")
summary(lm(y ~ x))

### **Ratio of non-synonymous and synonymous polymorphism, $\pi_N$/$\pi_S$**
#### Now we look at the correlation of the ratio $\pi_N$/$\pi_S$ with PC1
x <- pc1_mean
y <- polym$`Mean_piN/piS`[corr]
plot(x,y, type="n", xlab="PC1" , ylab=expression(italic(pi)[N]/italic(pi)[S]), xlim= c(min(x)-0.15, max(x)+0.15), ylim=c(0.05, 0.22), pch=18, col=mycol[match(species,aeg$species_code)], cex=1.4, cex.lab=1.2)
plotCI(x,y, li=polym$piNpiS_CI_lower[corr], ui=polym$piNpiS_CI_upper[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=19, cex=1, lwd=2)
# Add a weight to the regression (inverse of the variance of the estimator of piS)
# w <- 1 / (polym$piNpiS_CI_upper[corr] -polym$piNpiS_CI_lower[corr])^2
w <- NULL # unweighted regression
abline(lm(y ~ x,weights = w), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x,weights = w))$r.squared
summary(lm(y ~ x,weights = w))
text(x=x[-c(5,13)] - 0.35, y=y[-c(5,13)], labels=species[-c(5,13)], cex=1, col=mycol[match(species[-c(5,13)],aeg$species_code)]) # labels on the left
text(x=x[c(5,13)] + 0.35, y=y[c(5,13)], labels=species[c(5,13)], cex=1, col=mycol[match(species[c(5,13)],aeg$species_code)]) # labels on the right
legend(min(x)-0.5,y=max(y)+0.01, legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n", cex=1.4)
legend(min(x)-0.5,y=max(y)-0.01, legend = bquote("(With phylogenetic correction"~R^2 == "0.34*)"), bty = "n", cex=1) # with independent contrast correction

# The pic function is applied to the two vectors to be compared, then analyses can be done with the contrast values instead of the raw values.
x <- pc1_mean
y <- polym$`Mean_piN/piS`[corr]
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(y,tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
#w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
plot(x.pic,y.pic)
abline(lm(y.pic~x.pic-1,weights = w.pic))
summary(lm(y.pic~x.pic-1, weights = w.pic))

#### The relationships between $\pi_N$/$\pi_S$ and PC1 is positive and significant (Adjusted R-squared: 0.4364, p-value: 0.008331), even after correcting for the phylogeny (Adjusted R-squared:  0.289, p-value: 0.04141)
#### As for $\pi_S$, his relationship doesn't change if the median of PC1 per species is used instead of the mean. Not shown.
#The same correlation but with the median of PC1 per species
x <- pc1_median
y <- polym$`Mean_piN/piS`[corr]
plot(x,y, type="n", xlab="PC1", ylab=expression(italic(pi)[N]/italic(pi)[S]), xlim= c(min(x)-0.15, max(x)+0.15), ylim=c(0.05, 0.22), cex=1.4, pch=18, col=mycol[match(species,aeg$species_code)])
plotCI(x,y, li=polym$piNpiS_CI_lower[corr], ui=polym$piNpiS_CI_upper[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=19, cex=1, lwd=2)
abline(lm(y ~ x), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x))$r.squared
text(x=x[-5] - 0.35, y=y[-5], labels=species[-5], cex=1, col=mycol[match(species[-5],aeg$species_code)]) # labels on the left
text(x=x[5] + 0.35, y=y[5], labels=species[5], cex=1, col=mycol[match(species[5],aeg$species_code)]) # labels on the right
legend(min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n", cex=1.4)
summary(lm(y ~ x))

#### To summarize: we observe that synonymous diversity and selection against deleterious mutations significantly correlate with the selfing syndrome described by the first axis of the phenotypic PCA and this relationship is gradual. According to expectations, genome-wide diversity and the efficacy of selection against slightly deleterious mutations are lower in self-compatible species, and this effect is stronger for higher rates of self-fertilization.


## Figure export

pdf(file="figures/main/PiS_piNpiS_pc1_final_palette.pdf", height = 14, width=9, pointsize = 18)
pc1 <- pca_individuals$li[,1]
pc1_mean <- tapply(pc1,aeg$species_code, mean)
corr <-  match(species,row.names(polym))
layout(matrix(c(1,2),ncol=1)) # create space for two plots
# plot1. log10(piS) versus PC1 morphology
op <- par(mar = c(2,4.5,2,1))
x <- pc1_mean
y <- log10(polym$Mean_piS[corr])
plot(x,y, type="n", xlab="", ylab=expression(italic(pi)[S]~(log[10]~scale)), xlim= c(min(x)-0.15, max(x)+0.15), 
     ylim = c(min(y),max(y)-0.05*max(y)), pch=18, cex.lab= 1.4, col=mycol[match(species,aeg$species_code)], yaxt="n") 
plotCI(x,y, li=log10(polym$piS_CI_lower[corr]), ui=log10(polym$piS_CI_upper[corr]), add=T,
       col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, lwd=2)
axis(side=2, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05),cex.axis=1, padj=1)
# Add a weight to the regression (inverse of the variance of the estimator of log(piS): var(piS)/mean(piS)^2
# w <- polym$Mean_piS[corr]/ (polym$piS_CI_upper[corr] -polym$piS_CI_lower[corr])^2
w <- NULL # unweighted regression
abline(lm(y ~ x,weights = w), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x))$r.squared
# To get the R2 of the pic regression
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(y,tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
#w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
text(x=x[-c(2:5,11)] - 0.35, y=y[-c(2:5,11)], labels=species[-c(2:5,11)], cex=1.2,
     col=mycol[match(species[-c(2:5,11)],aeg$species_code)]) # labels on the left
text(x=x[c(2:5,11)], y=y[c(2:5,11)] + 0.07, labels=species[c(2:5,11)], cex=1.2,
     col=mycol[match(species[c(2:5,11)],aeg$species_code)]) # labels on top
legend(min(x)-0.5, y=max(y)-1.05, legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n", cex=1.4)
legend(min(x)-0.5, y=max(y)-1.2, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"***)"), bty = "n", cex=1) # with independent contrast correction
# plot2. piNpiS versus PC1 morphology
op <- par(mar = c(4,4.5,1,1))
x <- pc1_mean
y <- polym$`Mean_piN/piS`[corr]
plot(x,y, xlab="PC1", ylab=expression(italic(pi)[N]/italic(pi)[S]), xlim=c(min(x)-0.15, max(x)+0.15),
     ylim=c(0.08, 0.22), pch=18, cex=1.4, cex.lab= 1.4, col=mycol[match(species,aeg$species)])
plotCI(x,y, li=polym$piNpiS_CI_lower[corr], ui=polym$piNpiS_CI_upper[corr], add=T, 
       col=mycol[match(species,aeg$species_code)], pch=18, lwd=2)
# Add a weight to the regression (inverse of the variance of the estimator of log(piS): var(piS)/mean(piS)^2
# w <- polym$Mean_piS[corr]/ (polym$piS_CI_upper[corr] -polym$piS_CI_lower[corr])^2
w <- NULL # unweighted regression
abline(lm(y ~ x,weights = w), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x))$r.squared
# To get the R2 of the pic regression
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(y,tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
#w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
R2 <- summary(lm(y ~ x))$r.squared
text(x=x[-c(5,13)] - 0.35, y=y[-c(5,13)], labels=species[-c(5,13)], cex=1.2, col=mycol[match(species[-c(5,13)],aeg$species_code)]) # labels on the left
text(x=x[c(5,13)] + 0.35, y=y[c(5,13)], labels=species[c(5,13)], cex=1.2, col=mycol[match(species[c(5,13)],aeg$species_code)]) # labels on
legend(min(x)-0.5,y=max(y)+0.015, legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n", cex=1.4)
legend(min(x)-0.5,y=max(y), legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"*)"), bty = "n", cex=1) # with independent contrast correction
par(op)  # - reset to default
dev.off()



# 3) Verification of the correlation between polymorphism estimates ($\pi_N$/$\pi_S$ and $\pi_S$) and _F_, and between $\pi_N$/$\pi_S$ and $\pi_S$ ####

#### If _F_ mainly estimate the degree of self-fertilization, we should observe lower neutral diversity and efficacy of selection against deleterious mutations for species with higher _F_. Also, a negative correlation is expected between $\pi_N$/$\pi_S$ and diversity.
layout(matrix(c(1,2),ncol=2))
# piS
x <- polym$`Mean_Fit_W&C`[1:13]
y <- log10(polym$Mean_piS[1:13])
plot(x,y, xlab=expression(italic(F)), ylab=expression(italic(pi)[S]~(log[10]~scale)), xlim=c(min(x)-0.05, max(x)+0.05), ylim = c(min(y), max(y) + 0.1), yaxt="n", col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, cex.lab= 1.2)
axis(side=2, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05),cex.axis=1, padj=1)
abline(lm(y ~ x), lty=2 )
R2 <- summary(lm(y ~ x))$r.squared
text(x=x[-c(2:5,11:13)] - 0.05, y=y[-c(2:5,11:13)], labels=species[-c(2:5,11:13)], cex=1,
     col=mycol[match(species[-c(2:5,11:13)],aeg$species_code)]) # labels on the left
text(x=x[c(2:5,11)], y=y[c(2:5,11)] + 0.07, labels=species[c(2:5,11)], cex=1, 
     col=mycol[match(species[c(2:5,11)],aeg$species_code)]) # labels on top
text(x=x[c(12,13)], y=y[c(12,13)] - 0.07, labels=species[c(12,13)], cex=1, 
     col=mycol[match(species[c(12,13)],aeg$species_code)]) # labels on the bottom
legend(0.2,y=log10(0.003), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.4)
# piN/piS
x <- polym$`Mean_Fit_W&C`[1:13]
y <- polym$`Mean_piN/piS`[1:13]
plot(x,y, xlab=expression(italic(F)),ylab=expression(italic(pi)[N]/italic(pi)[S]), xlim= c(min(x)-0.05, max(x)+0.05), ylim = c(min(y), max(y) + 0.01),
     col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, cex.lab= 1.2) #, yaxt="n")
#axis(side=2, at=y, labels=y, las=2, cex.axis=0.8)
#axis(side = 2, at=y, labels = polym[na.omit(corr),2])
abline(lm(y ~ x), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x))$adj.r.squared
text(x=x, y=y+0.005, labels=species, cex=1, col=mycol[match(species,aeg$species_code)]) 
legend(min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex = 1.4)

#Figure export
pdf(file="figures/sup_mat/Suppl_Figure_Sx_Fis_piS_piNpiS.pdf", pointsize=14, height=9, width=18)
layout(matrix(c(1,2),ncol=2))
# piS
par(mar = c(4,4.5,2,1))
x <- polym$`Mean_Fit_W&C`[1:13]
y <- log10(polym$Mean_piS[1:13])
plot(x,y, xlab=expression(italic(F)), ylab=expression(italic(pi)[S]~(log[10]~scale)), xlim=c(min(x)-0.05, max(x)+0.05), ylim = c(min(y), max(y) + 0.1),
     yaxt="n", col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, cex.lab= 1.4)
axis(side=2, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05), cex.axis=1, padj=1)
abline(lm(y ~ x), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x))$r.squared
text(x=x[-c(2:5,11:13)] - 0.06, y=y[-c(2:5,11:13)], labels=species[-c(2:5,11:13)], cex=1,
     col=mycol[match(species[-c(2:5,11:13)],aeg$species_code)]) # labels on the left
text(x=x[c(2:5,11)], y=y[c(2:5,11)] + 0.07, labels=species[c(2:5,11)], cex=1, 
     col=mycol[match(species[c(2:5,11)],aeg$species_code)]) # labels on top
text(x=x[c(12,13)], y=y[c(12,13)] - 0.08, labels=species[c(12,13)], cex=1, 
     col=mycol[match(species[c(12,13)],aeg$species_code)]) # labels on the bottom
legend(0.2,y=log10(0.003), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.4)
# piN/piS
x <- polym$`Mean_Fit_W&C`[1:13]
y <- polym$`Mean_piN/piS`[1:13]
plot(x,y, xlab=expression(italic(F)),ylab=expression(italic(pi)[N]/italic(pi)[S]), xlim=c(min(x)-0.05, max(x)+0.05), ylim=c(0.09,0.22),
     col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, cex.lab= 1.4)
abline(lm(y ~ x), lty=2, lwd=1.5)
R2 <- summary(lm(y ~ x))$r.squared
text(x=x, y=y+0.005, labels=species, cex=0.8, col=mycol[match(species,aeg$species_code)]) 
legend(min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.4)
dev.off()

# piS
x <- log10(polym$Mean_piS[1:13])
y <- polym$`Mean_piN/piS`[1:13]
plot(x,y, xlab=expression(italic(pi)[S]~(log[10]~scale)),
     ylab=expression(italic(pi)[N]/italic(pi)[S]),
     xlim=c(min(x)-0.02, max(x)+0.02), ylim = c(min(y), max(y) + 0.01), 
     xaxt="n", col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, cex.lab= 1.4)
axis(side=1, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05), cex.axis=1.2, padj=0)
abline(lm(y ~ x), lty=2, lwd=1.5 )
R2 <- summary(lm(y ~ x))$r.squared
text(x=x, y=y+0.008, labels=species, cex=1, col=mycol[match(species,aeg$species_code)])
legend(-2,y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.4)
#### As expected, we observe that neutral diversity and the efficacy of selection against deleterious mutations are lower for species with higher _F_. Also, a negative correlation is observed between $\pi_N$/$\pi_S$ and diversity.

# Export figure
pdf(file="figures/sup_mat/Suppl_Figure_Sx_piS_piNpiS.pdf", pointsize=14, height=9, width=9)
x <- log10(polym$Mean_piS[1:13])
y <- polym$`Mean_piN/piS`[1:13]
plot(x,y, xlab=expression(italic(pi)[S]~(log[10]~scale)), ylab=expression(italic(pi)[S]/italic(pi)[N]), xlim=c(min(x)-0.02, max(x)+0.02), ylim = c(min(y), max(y) + 0.01), 
     xaxt="n", col=mycol[match(species,aeg$species_code)],pch=18, cex=1.4, cex.lab=1.4)
axis(side=1, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05), cex.axis=1, padj=0)
abline(lm(y ~ x), lty=2, lwd=1.5 )
R2 <- summary(lm(y ~ x))$r.squared
text(x=x, y=y+0.005, labels=species, cex=1, col=mycol[match(species,aeg$species_code)])
legend(log10(0.007),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.4)
dev.off()


# 4) Correlation with species range ####

layout(matrix(c(1,2),ncol=2))
op <- par(mar = c(5,4.5,2,1))
# piS
x <- polym$Species_range_grid_1.0[corr]
y <- polym$Mean_piS[corr]
# Add a weight
# w <- polym$Mean_piS[corr]^2 / (polym$piS_CI_upper[corr] - polym$piS_CI_lower[corr])^2
w <- NULL # unweighted regression
plot(x,y, xlab=expression("Species range (in 1000 km2)"), ylab=expression(italic(pi)[S]~(log[10]~scale)), pch=16, xlim=c(min(x)*0.95, max(x)*1.05), ylim = c(min(y), max(y)*1.2), col=mycol[match(species,aeg$species_code)], log="xy")
plotCI(x,y, li=log(polym$piS_CI_lower[corr]), ui=log(polym$piS_CI_upper[corr]), add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1.4, lwd=2)
coefreg <- lm(log(y) ~ log(x))$coefficients
ypred <- exp(coefreg[1])*c(1:max(x)*2)^coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2)
R2 <- summary(lm(log(y) ~ log(x),weights = w))$r.squared
text(x=x, y=1.2*y, labels=species, cex=1, col=mycol[match(species,aeg$species_code)])
legend(x=min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.4)
# piN/piS
x <- polym$Species_range_grid_1.0[corr]
y <- polym$`Mean_piN/piS`[corr]
# Add a weight
# w <- 1 / (polym$piNpiS_CI_upper[corr] - polym$piNpiS_CI_lower[corr])^2
w <- NULL # unweighted regression
plot(x,y, xlab=expression("Species range (in 1000 km2)"),ylab=expression(italic(pi)[N]/italic(pi)[S]), xlim= c(min(x)-0.05, max(x)+0.05), ylim = c(min(y), max(y) + 0.01), pch=18, cex=1.4, col=mycol[match(species,aeg$species_code)], log="x") #, yaxt="n")
plotCI(x,y, li=polym$piNpiS_CI_lower[corr], ui=polym$piNpiS_CI_upper[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
#axis(side=2, at=y, labels=y, las=2, cex.axis=0.8)
#axis(side = 2, at=y, labels = polym[na.omit(corr),2])
coefreg <- lm(y ~ log(x),weights = w)$coefficients
ypred <- coefreg[1]+log(c(1:max(x)*2))*coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2)
R2 <- summary(lm(y ~ log(x),weights = w))$r.squared
text(x=x[-9]*1.3, y=y[-9], labels=species[-9], cex=1, col=mycol[match(species[-9],aeg$species_code)])  # labels on teh right
text(x=x[9]-300, y=y[9], labels=species[9], cex=1, col=mycol[match(species[9],aeg$species_code)]) # Ata on the left
legend(min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n", cex=1.4)
par(op)

#Species range has no significant effect, and for $\pi_N/\pi_S$ the tendency is even opposite to the expectation.
#If we take mating systems into account, the effect of species range remains non-significant.

# Export figure
pdf(file="figures/main/piS_piNpiS_speciesrange_final_palette.pdf", height = 14, width=9, pointsize = 18)
layout(matrix(c(1,2),ncol=1))
# piS
op <- par(mar = c(2,4.5,2,1))
x <- polym$Species_range_grid_1.0[1:13]
y <- log10(polym$Mean_piS[1:13])
plot(x,y, yaxt="n",xlab=expression("Species range (in 1000 km2)"), log="x", ylab=expression(italic(pi)[S]~(log[10]~scale)), pch=18, xlim=c(min(x)*0.8, max(x)*1.05), ylim = c(min(y)-0.05, max(y)+0.05), col=mycol[match(species,aeg$species_code)], cex=1.4, cex.lab=1.4)
plotCI(x = x,y = y, li=log10(polym$piS_CI_lower[corr]), ui=log10(polym$piS_CI_upper[corr]), add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1.2, lwd=2)
coefreg <- lm(y ~ log(x))$coefficients
ypred <- coefreg[1]+log(c(1:max(x)*2))*coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2, lwd=1.5)
R2 <- summary(lm(y ~ log(x)))$r.squared
axis(side=2, at=log(c(0.001,0.002,0.005,0.01,0.02,0.05),10), labels=c(0.001,0.002,0.005,0.01,0.02,0.05), las=2, cex.axis=1)
# To get the R2 of the pic regression
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(y,tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
#w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
text(x=0.85*x, y=y+0.05, labels=species, cex=1.2, col=mycol[match(species,aeg$species_code)])
legend(x=0.7*min(x),y=max(y)+0.05, legend = bquote(R^2 == .(round(R2,2))), bty = "n",cex=1.2)
legend(x=0.7*min(x),y=max(y)-0.1, legend = bquote("(With phylogenetic correction "~R^2 == .(round(R2.pic,2))~")"), cex=1, bty = "n")
# piN/piS
op <- par(mar = c(4,4.5,1,1))
x <- polym$Species_range_grid_1.0[1:13]
y <- polym$`Mean_piN/piS`[1:13]
plot(x,y, xlab=expression("Species range (in 1000 km2)"),ylab=expression(italic(pi)[N]/italic(pi)[S]), xlim= c(0.8*min(x), 1.05*max(x)), ylim = c(0.9*min(y), max(y)*1.2), pch=18, cex=1.4, cex.lab=1.4,  col=mycol[match(species,aeg$species_code)], log="x") #, yaxt="n")
plotCI(x,y, li=polym$piNpiS_CI_lower[corr], ui=polym$piNpiS_CI_upper[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, lwd=2)
#axis(side=2, at=y, labels=y, las=2, cex.axis=0.8)
#axis(side = 2, at=y, labels = polym[na.omit(corr),2])
coefreg <- lm(y ~ log(x))$coefficients
ypred <- coefreg[1]+log(c(1:max(x)*2))*coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2, lwd=1.5)
R2 <- summary(lm(y ~ log(x)))$r.squared
# To get the R2 of the pic regression
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(y,tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
#w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
text(x=0.8*x, y=y*1.05, labels=species, cex=1.2, col=mycol[match(species,aeg$species_code)]) 
legend(0.7*min(x),y=max(y)*1.25, legend = bquote(R^2 == .(round(R2,2))), bty = "n",cex=1.4)
legend(0.7*min(x),y=max(y)*1.16, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~")"), bty = "n",cex=1)
par(op)  #- reset to default
dev.off()


# Same analysis but with the control of MS

layout(matrix(c(1,2),ncol=2))
pc1 <- pca_individuals$li[,1]
#sps <- as.character(aeg$species)
pc1_mean <- tapply(pc1,aeg$species, mean)
corr <-  match(species,row.names(polym))
x <- polym$Species_range_grid_1.0[1:13]
# y is the residual of the regression of piS on pc1
y <- lm(polym$Mean_piS[1:13]~pc1_mean[corr])$res
plot(x,y, xlab=expression("Species range (in 1000 km2)"), ylab=expression("Residual"), pch=18, xlim=c(min(x)*0.95, max(x)*1.05), ylim = c(min(y), max(y)*1.05), col=mycol[match(species,aeg$species)], log="x")
coefreg <- lm(y ~ log(x))$coefficients
ypred <- coefreg[1]+log(c(1:max(x)*2))*coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2)
R2 <- summary(lm(y ~ log(x)))$r.squared
text(x=x, y=y+0.0005, labels=species, cex=0.8, col=mycol[match(species,aeg$species)])
legend(x=min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n")
x <- polym$Species_range_grid_1.0[1:13]
# y is the residual of the regression of piN/piS on pc1
y <- lm(polym$`Mean_piN/piS`[1:13]~pc1_mean[corr])$res
plot(x,y, xlab=expression("Species range (in 1000 km2)"), ylab=expression("Residual"), pch=18, xlim=c(min(x)-0.05, max(x)*1.05), ylim = c(min(y), max(y)*1.05), col=mycol[match(species,aeg$species)], log="x")
coefreg <- lm(y ~ log(x))$coefficients
ypred <- coefreg[1]+log(c(1:max(x)*2))*coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2)
R2 <- summary(lm(y ~ log(x)))$r.squared
text(x=x, y=y+0.003, labels=species, cex=0.8, col=mycol[match(species,aeg$species)])
legend(x=min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n")
# Linear model with the two effects
mating <- pc1_mean[which(!is.na(corr))]
range <- polym$Species_range_grid_1.0[1:13]
piS <- log(polym$Mean_piS[1:13])
piNpiS <- polym$`Mean_piN/piS`[1:13]
summary(lm(piS~mating+range))
summary(lm(piNpiS~mating+range))

# Figure export
pdf(file="figures/sup_mat/Suppl_Figure_Sxx_piS_piNpiS_speciesrange_MScontrol_palette.pdf", height = 6, width=12, pointsize = 18)
layout(matrix(c(1,2),ncol=2))
op <- par(mar = c(4,4,2,2))
pc1 <- pca_individuals$li[,1]
#sps <- as.character(aeg$species)
pc1_mean <- tapply(pc1,aeg$species, mean)
corr <-  match(species,row.names(polym))
x <- polym$Species_range_grid_1.0[1:13]
# y is the residual of the regression of piS on pc1
y <- lm(polym$Mean_piS[1:13]~pc1_mean[corr])$res
plot(x,y, xlab=expression("Species range (in 1000 km2)"), ylab=expression("Residual"), pch=18, xlim=c(min(x)*0.95, max(x)*1.05), ylim = c(min(y), max(y)*1.05), col=mycol[match(species,aeg$species)], log="x")
coefreg <- lm(y ~ log(x))$coefficients
ypred <- coefreg[1]+log(c(1:max(x)*2))*coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2, lwd=1.5)
R2 <- summary(lm(y ~ log(x)))$r.squared
text(x=x, y=y+0.0005, labels=species, cex=0.8, col=mycol[match(species,aeg$species)])
legend(x=min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n")
x <- polym$Species_range_grid_1.0[1:13]
# y is the residual of the regression of piN/piS on pc1
y <- lm(polym$`Mean_piN/piS`[1:13]~pc1_mean[corr])$res
plot(x,y, xlab=expression("Species range (in 1000 km2)"), ylab=expression("Residual"), pch=18, xlim=c(min(x)*0.95, max(x)*1.05), ylim = c(min(y), max(y)*1.05), col=mycol[match(species,aeg$species)], log="x")
coefreg <- lm(y ~ log(x))$coefficients
ypred <- coefreg[1]+log(c(1:max(x)*2))*coefreg[2]
lines(c(1:max(x)*2),ypred,lty=2, lwd=1.5)
R2 <- summary(lm(y ~ log(x)))$r.squared
text(x=x, y=y+0.003, labels=species, cex=0.8, col=mycol[match(species,aeg$species)])
legend(x=min(x),y=max(y), legend = bquote(R^2 == .(round(R2,2))), bty = "n")
dev.off()



# 5) Correlation with linked selection parameters ####

#### Plot of the maximum _$\pi[S]$_ obtained from the fit of the lonked selection model

pc1 <- pca_individuals$li[,1]
#sps <- as.character(aeg$species)
pc1_mean <- tapply(pc1,aeg$species, mean)
corr <-  match(species,row.names(polym))
#layout(matrix(c(1,2,3,4),ncol=2))
# pi_max/piS
x <- pc1_mean
y <- polym$piMax_fitted[corr]
plot(x,y,
     xlab="PC1",
     ylab=expression(italic(pi)[max]),
     xlim=c(min(x)-0.15, max(x)+0.15),
     ylim=c(min(y)*0.3,max(y)*1.5),
     pch=18, cex=1.4,
     log="y",
     col=mycol[match(species,aeg$species)]) #, yaxt="n")
plotCI(x,y, li=polym$piMax_fitted_inf[corr], ui=polym$piMax_fitted_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
abline(lm(log10(y) ~ x), lty=2, lwd=1.5)
R2 <- summary(lm(log10(y) ~ x))$r.squared
# To get the R2 of the pic regression
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(y,tree)

# Figure export
pdf(file="figures/main/LinkedSelection_pc1_final_palette.pdf", height = 6, width = 14, pointsize = 18)
layout(matrix(c(1,2),ncol=2))
pc1 <- pca_individuals$li[,1]
#sps <- as.character(aeg$species)
pc1_mean <- tapply(pc1,aeg$species, mean)
corr <-  match(species,row.names(polym))
#layout(matrix(c(1,2,3,4),ncol=2))
# pi_max/piS
x <- pc1_mean
y <- polym$piMax_fitted[corr]
plot(x,y,
     xlab="PC1",
     ylab=expression(italic(pi)[max]),
     xlim=c(min(x)-0.15, max(x)+0.15),
     ylim=c(min(y)*0.3,max(y)*1.5),
     pch=18, cex=1.4,cex.lab=1.2,
     log="y",
     col=mycol[match(species,aeg$species)]) #, yaxt="n")
plotCI(x,y, li=polym$piMax_fitted_inf[corr], ui=polym$piMax_fitted_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
abline(lm(log10(y) ~ x), lty=2, lwd=1.5)
R2 <- summary(lm(log10(y) ~ x))$r.squared
# To get the R2 of the pic regression
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(x,tree)
y.pic <- pic(log(y),tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
#w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
#text(x=x-0.3, y=y+0.05, labels=species, cex=0.8,col=mycol[match(species,aeg$species)]) 
text(x=x[-c(5,12,13)] - 0.35, y=y[-c(5,12,13)], labels=species[-c(5,12,13)], cex=0.8, col=mycol[match(species[-c(5,12,13)],aeg$species_code)]) # labels on the left
text(x=x[c(5,12,13)] + 0.35, y=y[c(5,12,13)], labels=species[c(5,12,13)], cex=0.8, col=mycol[match(species[c(5,12,13)],aeg$species_code)]) # labels on the right
legend(x=min(x)-0.6,y=3, legend = bquote(R^2 == .(round(R2,3))), bty = "n",cex=1)
legend(min(x)-0.6,y=1.5, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,3))~")"), cex=0.8,bty = "n")

corr <-  match(species,row.names(polym))
#layout(matrix(c(1,2,3,4),ncol=2))
# pi_max/piS
x <- polym$Species_range_grid_1.0[corr]
y <- polym$piMax_fitted[corr]
plot(x,y,
     xlab="Species range (in 1000 km2)",
     ylab=expression(italic(pi)[max]),
     xlim=c(min(x)-0.15, max(x)+0.15),
     ylim=c(min(y)*0.3,max(y)*1.5),
     pch=18, cex=1.4, cex.lab=1.2,
     log="xy",
     col=mycol[match(species,aeg$species)]) #, yaxt="n")
plotCI(x,y, li=polym$piMax_fitted_inf[corr], ui=polym$piMax_fitted_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
abline(lm(log10(y) ~ log10(x)), lty=2, lwd=1.5)
R2 <- summary(lm(log10(y) ~ log10(x)))$r.squared
# To get the R2 of the pic regression
names(x) <- sort(newlabels)
names(y) <- sort(newlabels)
x.pic <- pic(log(x),tree)
y.pic <- pic(log(y),tree)
# The variance of contrasts can be obtained from the pic function. The variance depends of branch length
#w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
w.pic <- NULL # unweighted regression
# Note that the regression with contrast must be forced to pass through the origin
R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
#text(x=x-0.3, y=y+0.05, labels=species, cex=0.8,col=mycol[match(species,aeg$species)]) 
text(x=x[-c(5,12,13)] - 0.35, y=y[-c(5,12,13)], labels=species[-c(5,12,13)], cex=0.8, col=mycol[match(species[-c(5,12,13)],aeg$species_code)]) # labels on the left
text(x=x[c(5,12,13)] + 0.35, y=y[c(5,12,13)], labels=species[c(5,12,13)], cex=0.8, col=mycol[match(species[c(5,12,13)],aeg$species_code)]) # labels on the right
legend(x=150,y=5, legend = bquote(R^2 == .(round(R2,2))~"**"), bty = "n",cex=1)
legend(x=150,y=2.5, legend = bquote("(With phylogenetic"), cex=0.8,bty = "n")
legend(x=150,y=1.5, legend = bquote("correction"~R^2 == .(round(R2.pic,2))~")"), cex=0.8,bty = "n")

dev.off()



# #### After fitting a logistic curve on the relationship between polymorphism and recombination, we can extract two useful statistics: the slope of the relationship between polymorphism and recombination at the inflexion point (the lowest the slope, the strongest the linked selection effects), and the maximum polymorphism that could be reached for free recombination.
# 
# #### **Synonymous polymorphism maximum _$\pi[S]$_**
# pc1 <- pca_individuals$li[,1]
# #sps <- as.character(aeg$species)
# pc1_mean <- tapply(pc1,aeg$species, mean)
# corr <-  match(species,row.names(polym))
# layout(matrix(c(1,2,3,4),ncol=2))
# # pi_max/piS
# x <- pc1_mean
# y <- polym$Ratio_piS_piMax[corr]
# plot(x,y, xlab="PC1",ylab=expression(italic(pi)[max]*"/"*italic(pi)[S]), xlim= c(min(x)-0.15, max(x)+0.15), ylim = c(1,6.5), pch=18, cex=1.4, col=mycol[match(species,aeg$species)]) #, yaxt="n")
# plotCI(x,y, li=polym$Ratio_piS_piMax_inf[corr], ui=polym$Ratio_piS_piMax_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
# abline(lm(y ~ x), lty=2, lwd=1.5)
# R2 <- summary(lm(y ~ x))$adj.r.squared
# # To get the R2 of the pic regression
# names(x) <- sort(newlabels)
# names(y) <- sort(newlabels)
# x.pic <- pic(x,tree)
# y.pic <- pic(y,tree)
# # The variance of contrasts can be obtained from the pic function. The variance depends of branch length
# #w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
# w.pic <- NULL # unweighted regression
# # Note that the regression with contrast must be forced to pass through the origin
# R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
# text(x=x-0.3, y=y+0.05, labels=species, cex=0.8,col=mycol[match(species,aeg$species)]) 
# legend(x=min(x)-0.5,y=6.5, legend = bquote(R^2 == .(round(R2,2))~"**"), bty = "n",cex=1.4)
# legend(min(x)-0.5,y=5.7, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"**)"), bty = "n")
# # Slope
# x <- pc1_mean
# y <- polym$Slope_rec_piS[corr]
# plot(x,y, xlab="PC1",ylab="Slope", xlim= c(min(x)-0.15, max(x)+0.15), ylim = c(0, 0.06),
#      pch=18, cex=1.4, col=mycol[match(species,aeg$species)]) #, yaxt="n")
# plotCI(x,y, li=polym$Slope_rec_piS_inf[corr], ui=polym$Slope_rec_piS_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
# abline(lm(y ~ x), lty=2, lwd=1.5)
# R2 <- summary(lm(y ~ x))$adj.r.squared
# # To get the R2 of the pic regression
# names(x) <- sort(newlabels)
# names(y) <- sort(newlabels)
# x.pic <- pic(x,tree)
# y.pic <- pic(y,tree)
# # The variance of contrasts can be obtained from the pic function. The variance depends of branch length
# #w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
# w.pic <- NULL # unweighted regression
# # Note that the regression with contrast must be forced to pass through the origin
# R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
# text(x=x+0.25, y=y+0.002, labels=species, cex=0.8,col=mycol[match(species,aeg$species)]) 
# legend(x=-4,y=0.06, legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n",cex=1.4)
# legend(x=-4,y=0.0525, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"***)"), bty = "n")
# # pi_max
# x <- pc1_mean
# y <- polym$Max_piS[corr]
# plot(x,y, xlab="PC1",ylab=expression(italic(pi)[max]*"/"*italic(pi)[S]), xlim= c(min(x)-0.15, max(x)+0.15), ylim = c(0,0.05), pch=18, cex=1.4, col=mycol[match(species,aeg$species)]) #, yaxt="n")
# plotCI(x,y, li=polym$Max_piS_inf[corr], ui=polym$Max_piS_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
# abline(lm(y ~ x), lty=2, lwd=1.5)
# R2 <- summary(lm(y ~ x))$adj.r.squared
# # To get the R2 of the pic regression
# names(x) <- sort(newlabels)
# names(y) <- sort(newlabels)
# x.pic <- pic(x,tree)
# y.pic <- pic(y,tree)
# # The variance of contrasts can be obtained from the pic function. The variance depends of branch length
# #w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
# w.pic <- NULL # unweighted regression
# # Note that the regression with contrast must be forced to pass through the origin
# R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
# text(x=x-0.3, y=y+0.05, labels=species, cex=0.8,col=mycol[match(species,aeg$species)]) 
# legend(x=min(x)-0.5,y=0.05, legend = bquote(R^2 == .(round(R2,2))~"**"), bty = "n",cex=1.4)
# legend(min(x)-0.5,y=0.043, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"**)"), bty = "n")
# 
# # Figure export
# pdf(file="figures/main/LinkedSelection_pc1_final_palette.pdf", height = 15, width= 8, pointsize = 18)
# pc1 <- pca_individuals$li[,1]
# #sps <- as.character(aeg$species)
# pc1_mean <- tapply(pc1,aeg$species, mean)
# corr <-  match(species,row.names(polym))
# layout(matrix(c(1,2),ncol=1))
# # pi_max/piS
# op <- par(mar = c(2,4.5,2,1))
# x <- pc1_mean
# y <- polym$Ratio_piS_piMax[corr]
# plot(x,y, xlab="PC1",ylab=expression(italic(pi)[max]*"/"*italic(pi)[S]), xlim= c(min(x)-0.15, max(x)+0.15), ylim = c(1,6), pch=18, cex=1.4, cex.lab=1.4, col=mycol[match(species,aeg$species)]) #, yaxt="n")
# plotCI(x,y, li=polym$Ratio_piS_piMax_inf[corr], ui=polym$Ratio_piS_piMax_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
# abline(lm(y ~ x), lty=2, lwd=1.5)
# R2 <- summary(lm(y ~ x))$adj.r.squared
# # To get the R2 of the pic regression
# names(x) <- sort(newlabels)
# names(y) <- sort(newlabels)
# x.pic <- pic(x,tree)
# y.pic <- pic(y,tree)
# # The variance of contrasts can be obtained from the pic function. The variance depends of branch length
# #w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
# w.pic <- NULL # unweighted regression
# # Note that the regression with contrast must be forced to pass through the origin
# R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
# #text(x=x-0.3, y=y+0.05, labels=species, cex=0.8,col=mycol[match(species,aeg$species)]) 
# text(x=x[-c(5,12,13)] - 0.35, y=y[-c(5,12,13)], labels=species[-c(5,12,13)], cex=1.2, col=mycol[match(species[-c(5,12,13)],aeg$species_code)]) # labels on the left
# text(x=x[c(5,12,13)] + 0.35, y=y[c(5,12,13)], labels=species[c(5,12,13)], cex=1.2, col=mycol[match(species[c(5,12,13)],aeg$species_code)]) # labels on the right
# legend(x=min(x)-0.5,y=6.3, legend = bquote(R^2 == .(round(R2,2))~"**"), bty = "n",cex=1.4)
# legend(min(x)-0.5,y=5.7, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"**)"), bty = "n")
# # Slope
# op <- par(mar = c(4,4.5,2,1))
# x <- pc1_mean
# y <- polym$Slope_rec_piS[corr]
# plot(x,y, xlab="PC1",ylab="Slope", xlim= c(min(x)-0.15, max(x)+0.15), ylim = c(0, 0.045),
#      pch=18, cex=1.4, cex.lab=1.4, col=mycol[match(species,aeg$species)]) #, yaxt="n")
# plotCI(x,y, li=polym$Slope_rec_piS_inf[corr], ui=polym$Slope_rec_piS_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
# abline(lm(y ~ x), lty=2, lwd=1.5)
# R2 <- summary(lm(y ~ x))$adj.r.squared
# # To get the R2 of the pic regression
# names(x) <- sort(newlabels)
# names(y) <- sort(newlabels)
# x.pic <- pic(x,tree)
# y.pic <- pic(y,tree)
# # The variance of contrasts can be obtained from the pic function. The variance depends of branch length
# #w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
# w.pic <- NULL # unweighted regression
# # Note that the regression with contrast must be forced to pass through the origin
# R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
# text(x=x[-c(1,2,6,9:11,13)]+0.35, y=y[-c(1,2,6,9:11,13)]+0.002, labels=species[-c(1,2,6,9:11,13)], cex=1.2,col=mycol[match(species[-c(1,2,6,9:11,13)],aeg$species)]) 
# text(x=x[c(1,2,6,10)] - 0.4, y=y[c(1,2,6,10)], labels=species[c(1,2,6,10)], cex=1.2, col=mycol[match(species[c(1,2,6,10)],aeg$species_code)]) # labels on the left
# text(x=x[c(13)], y=y[c(13)]- 0.002, labels=species[c(13)], cex=1.2, col=mycol[match(species[c(13)],aeg$species_code)]) # Tur below
# text(x=x[c(11)], y=y[c(11)]+ 0.002, labels=species[c(11)], cex=1.2, col=mycol[match(species[c(11)],aeg$species_code)]) # Aun above
# text(x=x[c(9)]-0.2, y=y[c(9)]+ 0.002, labels=species[c(9)], cex=1.2, col=mycol[match(species[c(9)],aeg$species_code)]) # Ata, left-above
# legend(x=-4.5,y=0.008, legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n",cex=1.4)
# legend(x=-4.5,y=0.003, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"***)"), bty = "n")
# #text(x=x[c(2:5,11)], y=y[c(2:5,11)] + 0.07, labels=species[c(2:5,11)], cex=1.2,
# #     col=mycol[match(species[c(2:5,11)],aeg$species_code)]) # labels on top
# #legend(min(x)-0.5, y=max(y)-1.05, legend = bquote(R^2 == .(round(R2,2))~"***"), bty = "n", cex=1.4)
# #legend(min(x)-0.5, y=max(y)-1.2, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"***)"), bty = "n", cex=1) # with independent contrast correction
# dev.off()
# 
# pdf(file="figures/sup_mat/Suppl_Figure_Sxx_pimax_vs_pc1.pdf", height = 8, width= 12, pointsize = 18)
# op <- par(mar = c(4,4.5,2,1))
# x <- pc1_mean
# y <- polym$Max_piS[corr]
# plot(x,y, xlab="PC1",ylab=expression(italic(pi)[max]), xlim= c(min(x)-0.15, max(x)+0.15), ylim = c(0,0.05), pch=18, cex=1.4, col=mycol[match(species,aeg$species)]) #, yaxt="n")
# plotCI(x,y, li=polym$Max_piS_inf[corr], ui=polym$Max_piS_sup[corr], add=T, col=mycol[match(species,aeg$species_code)], pch=18, cex=1, lwd=2)
# abline(lm(y ~ x), lty=2, lwd=1.5)
# R2 <- summary(lm(y ~ x))$adj.r.squared
# # To get the R2 of the pic regression
# names(x) <- sort(newlabels)
# names(y) <- sort(newlabels)
# x.pic <- pic(x,tree)
# y.pic <- pic(y,tree)
# # The variance of contrasts can be obtained from the pic function. The variance depends of branch length
# #w.pic <- 1/pic(y,tree,var.contrasts = T)[,2]
# w.pic <- NULL # unweighted regression
# # Note that the regression with contrast must be forced to pass through the origin
# R2.pic <- summary(lm(y.pic~x.pic-1, weights = w.pic))$r.squared
# text(x=x-0.3, y=y+0.05, labels=species, cex=0.8,col=mycol[match(species,aeg$species)]) 
# legend(x=min(x)-0.5,y=0.05, legend = bquote(R^2 == .(round(R2,2))~"**"), bty = "n",cex=1.4)
# legend(min(x)-0.5,y=0.043, legend = bquote("(With phylogenetic correction"~R^2 == .(round(R2.pic,2))~"**)"), bty = "n")
# dev.off()
