# Analysis of morphofunctional traits in Aegilops/Triticum species
# Part 2. Trait imputation and distribution
# author: Concetta Burgarella (Uppsala University)
# date: September 2021, updated November 2022


#packages to install
if (!require("psych"))   install.packages("psych", dependencies = TRUE)
if (!require("missMDA"))   install.packages("missMDA", dependencies = TRUE)
if (!require("FactoMineR"))   install.packages("FactoMineR", dependencies = TRUE)
if (!require("dplyr"))   install.packages("dplyr", dependencies = TRUE)

# 1. Diversity and distribution of traits before imputation #### 


# load the raw (non imputed) measures
non_imputed_morpho <- read.csv("outputs/sup_mat/Supplementary_table_Sx_Morpho_traits_original.csv")

# This table includes the measures directly taken on the plants. 
# List of species
non_imputed_morpho$species

# List the variables 
colnames(non_imputed_morpho)


#### **Boxplots of each trait per species**
# Create vector of colors for the 13 Aegilops/Triticum species
mycol <- rep(NA, nrow(non_imputed_morpho))
mycol[c(which(non_imputed_morpho$species_code == "Amu"))] <- "#085CF8"
mycol[c(which(non_imputed_morpho$species_code == "Asp"))] <- "#1684A6"   # "royalblue"
mycol[c(which(non_imputed_morpho$species_code == "Ash"))] <- "#3C9E49"
mycol[c(which(non_imputed_morpho$species_code == "Aca"))] <- "#65AF1E"
mycol[c(which(non_imputed_morpho$species_code == "Alo"))] <- "#98BB18"
mycol[c(which(non_imputed_morpho$species_code == "Aum"))] <- "#C7C612" 
mycol[c(which(non_imputed_morpho$species_code == "Abi"))] <- "#F3CC1D"
mycol[c(which(non_imputed_morpho$species_code == "Aco"))] <- "#FDAF55"
mycol[c(which(non_imputed_morpho$species_code == "Tbo"))] <- "#FE8F7B"
mycol[c(which(non_imputed_morpho$species_code == "Ase"))] <- "#FC6A9B"
mycol[c(which(non_imputed_morpho$species_code == "Aun"))] <- "#F64497" 
mycol[c(which(non_imputed_morpho$species_code == "Ata"))] <- "#E92952"
mycol[c(which(non_imputed_morpho$species_code == "Tur"))] <- "#D70500"



# 2. Imputation of missing data ####

#### For ovary and stigma only the first 3 measures per individual were included because the 4th is available only for few individuals.
toimpute <- non_imputed_morpho[, -c(1:3, 25, 29)]       # exclude columns
# Traits to be imputed are:
colnames(toimpute)

sum(is.na(toimpute)) # missing values over 2862 cells
sum(is.na(toimpute)) / (dim(toimpute)[1]*dim(toimpute)[2]) # percentage of missing values 

# How many samples have complete traits?
nrow(na.omit(toimpute)) # 39 samples are complete

#### Missing values are predicted using a principal component analysis. The optimum number of components to retain is estimated automatically as the number of PCs that minimizes the "mean squared error of prediction". \
# Estimation of the optimum nb of PC for imputation
nb <- estim_ncpPCA(toimpute,ncp.min=0,ncp.max=10, method.cv="Kfold")  
nb
# Plot the optimum nb of PC for imputation
plot(nb$criterion~names(nb$criterion),xlab="number of dimensions", ylab="mean square error of prediction")
#fig.cap="mean square error of prediction (MSEP) obtained for each number of PC by cross-validation approach"
# Perform imputation based on estimated npc
ncp <- nb$ncp
res.imp <- imputePCA(toimpute, ncp = ncp)    # single imputation
# Save as R object
save(list=ls()[which(ls()=="res.imp")], file = paste0("outputs/morphology/ImputePCA_results.Rdata"))

#Three PCs minimize the square error of prediction
#### Check that imputed values are positive values


#The imputation step can be peformed only once (for this set the eval option to TRUE in the previous chunk). Later, just load the imputed dataset.
# Load the imputed dataset
load("outputs/morphology/ImputePCA_results.Rdata")

# Check imputation
s <- summary(res.imp$completeObs) 
s

#### Some variables have negative values. We substitute negative values with 0 (since negative values have no biological meaning).
imp_obs_data <- cbind(non_imputed_morpho[,1:3], res.imp$completeObs)
# Substitute negative values with 0
imp_obs_data[,-(1:3)][imp_obs_data[,-(1:3)] < 0] <- 0

#### We now compare non imputed and imputed data by performing a PCA (in non imputed data, missing values are imputed by the mean of the variable).
# pca comparison between non imputed and imputed data
res.PCA.N <- PCA(toimpute, graph=F)
res.PCA <- PCA(res.imp$completeObs, graph=F)

# Imputed vs nonimputed data
plot(res.PCA.N,choix="var", title="non imputed measures",autoLab = "yes", cex=0.8)  # graph of the variables
plot(res.PCA.N,choix="ind", title="non imputed",autoLab = "yes")  # graph of the individuals
plot(res.PCA,choix="var",title="imputed measures", autoLab = "yes", cex=0.8)    # graph of the variables
plot(res.PCA,choix="ind",title="imputed", autoLab = "yes")    # graph of the individuals
#### With imputation, the relative position of some variable change and imputed variables have stronger contribution to axes. 


# 3. Calculation of useful summary statistics on imputed measures ####

#### We calculate summary statistics on imputed measures and use them as additional traits to characterize species mating system. 

# Statistic are: \
# * Eta_len_norm: mean anther length normalised by the flower length. \
# * Eta_wid_norm: mean anther width normalised by the flower length. \
# * Sti_norm: mean stigma length normalised by the flower length. \
# * Maleinvest: male investiment, calculated as mean anther length/mean ovary length. \
# * Compactness: a measure of spikelet overlap calculated as the ratio of (mean flower length* number_fertile_flower)/mean spikelet length. \
# * AF_index: autonomous seed set, a measure of the auto-fertilization rate, calculated as self-fertilised_seed_number/(self-fertilized_spikelet_number*number_fertile_flower/spikelet). \


# Calculate summary measures from imputed data
mean_eta_len <- apply(imp_obs_data[,10:15],MARGIN = 1, mean)     # mean of anther length
mean_eta_wid <- apply(imp_obs_data[,16:21],MARGIN = 1, mean)     # mean of anther width
eta_len_norm <- mean_eta_len/imp_obs_data$flower_length_mm  # mean anther lenght normalised by flower length
eta_wid_norm <- mean_eta_wid/imp_obs_data$flower_length_mm  # mean anther width normalised by flower length
sti_mean <- apply(imp_obs_data[,25:27],MARGIN = 1, mean)  # mean of stigma length
sti_norm <- sti_mean/imp_obs_data$flower_length_mm # mean of stigma length normalised by flower length
# maleinvest = taille des anthères normalisée par la taille des ovaires
ova_mean <- apply(imp_obs_data[,22:24],MARGIN = 1, mean)  # mean of ovario length
maleinvest <- mean_eta_len/ova_mean
# spikelet compactness = flower length/spikelet_length 
l <- imp_obs_data$spikelet_length_mm
compactness <- imp_obs_data$flower_length_mm/l
# autonomous seed set (following Escobar et al. 2010)
fertile_flo_x_spikelet_FR <- imp_obs_data$tot3_fertile/3
AF_index <- imp_obs_data$tot_AF_fertile/(imp_obs_data$n_AF_splets*fertile_flo_x_spikelet_FR)

# Check the values of the autonomous seed set 
summary(AF_index)

#### There are Inf and NaN values. Inf values are due to fertile_flo_x_spikelet_FR = 0, NaN values are due to both tot_AF_fertile and fertile_flo_x_spikelet_FR = 0. \
#### We substitute zero values of fertile_flo_x_spikelet_FR to 1 and recalculate the AF_index. This will make AF_index > 0 when tot_AF_fertile > 0, and AF_index = 0 when tot_AF_fertile = 0.


# Correct Inf and NaN values
fertile_flo_x_spikelet_FR[fertile_flo_x_spikelet_FR == 0] <- 1
AF_index <- imp_obs_data$tot_AF_fertile/(imp_obs_data$n_AF_splets*fertile_flo_x_spikelet_FR)
summary(AF_index)

#### We add the newly statistics to the table of imputed measures.
imputed_data <- cbind(imp_obs_data, mean_eta_len, mean_eta_wid, eta_len_norm, eta_wid_norm, sti_mean, ova_mean, sti_norm, maleinvest, compactness, AF_index)

# check imputed values
summary(imputed_data) # no negative values anymore


#### Export the table of imputed data. 
write.csv(imputed_data,"outputs/sup_mat/Supplementary_table_Sx_Morpho_traits_imputed.csv", quote = F, row.names = F)
save(list=ls()[which(ls()=="imputed_data")],file = "outputs/morphology/imputed_morpho_data.Rdata")


# 4. Plots of diversity and distribution of the imputed traits and summary statistics/traits build on them ####


# for anthers, stigmas and ovaries, include only mean values.

pdf(file = "figures/sup_mat/Suppl_Figure_Sx_boxplots_imputed_traits.pdf", pointsize = 8, height=10, width=6)
par(mfrow = c(6, 3))
par(cex = 0.6)
par(mar = c(2.2, 2.2, 2.5, 2.2), oma = c(1.2, 1, 0, 1))
for(i in c(4:9,28:ncol(imputed_data))){  
  boxplot(imputed_data[,i] ~ as.factor(imputed_data$species_code), 
          col=unique(mycol), main=colnames(imputed_data)[i], xlab="", ylab="", 
          las=2, cex.axis=1.1, cex.main=1.5)
}
dev.off()




#### **Pairwise correlation among all traits**

source("http://www.sthda.com/upload/rquery_cormat.r")
mydata <- imputed_data[,-c(1:3,10:27)]
require("corrplot")
rquery.cormat(mydata)

#### We look more in detail at the pairwise correlation of the summary statistics:

var <- c("eta_len_norm", "eta_wid_norm", "sti_norm", "ova_norm", "compactness", "AF_index", "maleinvest")
mydata2 <- imputed_data[,which(names(imputed_data) %in% var)]

library(PerformanceAnalytics)
chart.Correlation(mydata2, histogram = TRUE, method = "pearson")

#  export correlograms as pdf figure
pdf(file = paste("figures/sup_mat/Suppl_Figure_Sx_correlogram_imputed_traits.pdf",sep=""), pointsize = 25, height=10, width=15)
# plot1
mydata <- imputed_data[,-c(1:3,10:27,30)]
require("corrplot")
rquery.cormat(mydata)
dev.off()

pdf(file = paste("figures/sup_mat/Suppl_Figure_Sx_correlogram3_imputed_traits.pdf",sep=""), pointsize = 20, height=10, width=15)
# plot2
chart.Correlation(mydata2, histogram = TRUE, method = "pearson")
dev.off()

