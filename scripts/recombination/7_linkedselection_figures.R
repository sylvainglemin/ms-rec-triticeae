# Script to retrieve results of the linked selection model and combine them into figure
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022

if (!require("ggplot2"))   install.packages("ggplot2", dependencies = TRUE)

species_list <- c(
  "Ae_bicornis",
  "Ae_caudata",
  "Ae_comosa",
  "Ae_longissima",
  "Ae_mutica",
  "Ae_searsii",
  "Ae_sharonensis",
  "Ae_speltoides",
  "Ae_tauschii",
  "Ae_umbellulata",
  "Ae_uniaristata",
  "T_monococcum",
  "T_urartu"
)

spcode <- rep(NA, length(species_list))
spcode[c(which(species_list == "Ae_mutica"))] <- "Amu"
spcode[c(which(species_list == "Ae_speltoides"))] <- "Asp"   # "royalblue"
spcode[c(which(species_list == "Ae_sharonensis"))] <- "Ash"
spcode[c(which(species_list == "Ae_caudata"))] <- "Aca"
spcode[c(which(species_list == "Ae_longissima"))] <- "Alo"
spcode[c(which(species_list == "Ae_umbellulata"))] <- "Aum" 
spcode[c(which(species_list == "Ae_bicornis"))] <- "Abi"
spcode[c(which(species_list == "Ae_comosa"))] <- "Aco"
spcode[c(which(species_list == "T_monococcum"))] <- "Tmo"
spcode[c(which(species_list == "Ae_searsii"))] <- "Ase"
spcode[c(which(species_list == "Ae_uniaristata"))] <- "Aun" 
spcode[c(which(species_list == "Ae_tauschii"))] <- "Ata"
spcode[c(which(species_list == "T_urartu"))] <- "Tur"


mycol <- rep(NA, length(species_list))
mycol[c(which(species_list == "Ae_mutica"))] <- "#085CF8"
mycol[c(which(species_list == "Ae_speltoides"))] <- "#1684A6"   # "royalblue"
mycol[c(which(species_list == "Ae_sharonensis"))] <- "#3C9E49"
mycol[c(which(species_list == "Ae_caudata"))] <- "#65AF1E"
mycol[c(which(species_list == "Ae_longissima"))] <- "#98BB18"
mycol[c(which(species_list == "Ae_umbellulata"))] <- "#C7C612" 
mycol[c(which(species_list == "Ae_bicornis"))] <- "#F3CC1D"
mycol[c(which(species_list == "Ae_comosa"))] <- "#FDAF55"
mycol[c(which(species_list == "T_monococcum"))] <- "#FE8F7B"
mycol[c(which(species_list == "Ae_searsii"))] <- "#FC6A9B"
mycol[c(which(species_list == "Ae_uniaristata"))] <- "#F64497" 
mycol[c(which(species_list == "Ae_tauschii"))] <- "#E92952"
mycol[c(which(species_list == "T_urartu"))] <- "#D70500"



load("outputs/morphology/PCA_morphological_traits.Rdata")
pc1 <- pca_individuals$li[,1]
pc1_mean <- tapply(pc1,aeg$species_code, mean)

mydata <- c()
for(FILE in list.files(path = "outputs/recombination/",pattern = "Results_BSfit_filter300") ){
  mydata <-  rbind(mydata, read.table(paste("outputs/recombination/",FILE,sep=""),header=T))
}
mydata <- mydata[-6]
for(FILE in list.files(path = "outputs/recombination/",pattern = "Results_BSfit_FixedFis_filter300") ){
  mydata <-  rbind(mydata, read.table(paste("outputs/recombination/",FILE,sep=""),header=T))
}
# Replace T_boeticum by T_monococcum so also change species list
mydata$Species <- ifelse(mydata$Species=="T_boeticum","T_monococcum",mydata$Species)

mydata$spcode <- spcode[match(mydata$Species,species_list)]
mydata$pc1 <- pc1_mean[match(mydata$Species,species_list)]

mydata$ratio <- mydata$piS_est/mydata$piS_obs


ord <- order(aggregate(mydata$ratio,by = list(mydata$Species),median)$x)
mydata$Species <- factor(mydata$Species,levels = species_list[ord])



pdf("figures/sup_mat/linked-selection_model.pdf",width = 8,height = 6)
par(mar=c(10,6,4,4))
boxplot(mydata$ratio~mydata$Species,log="y",col=mycol[match(levels(mydata$Species),species_list)],xlab=NULL,ylab=NULL,
        main=expression(pi[max]/pi[S]),las=2)
dev.off()

mydata2 <- mydata[mydata$spcode!="Abi" & mydata$spcode!="Ash" & mydata$spcode!="Alo" & mydata$spcode!="Ase",]
mydatat2 <- droplevels(mydata2)
mydata3 <- aggregate(list(mydata2[,c("pc1","ratio")]),by=list(mydata2$Species),median)
mydata3$spcode <- aggregate(list(mydata2[,c("spcode")]),by=list(mydata2$spcode),function(x) as.character(x[1]))[,2]
names(mydata3) <- c("Species","pc1","ratio","spcode")

pdf("figures/sup_mat/linked-selection_model_goodspecies.pdf",width = 6,height = 4)
par(mar=c(5,5,3,4))
plot(mydata3$ratio~mydata3$pc1,
     ylim=c(1,max(mydata3$ratio)),
     xlab="PC1",ylab=expression(pi[max]/pi[S]),las=2,log="y")
points(mydata3$ratio~mydata3$pc1,pch=16,cex=0.8,col=mycol[match(mydata3$Species,species_list)],
     xlab="PC1",ylab=expression(pi[max]/pi[S]),las=2)
text(x=mydata3$pc1,y=mydata3$ratio*0.85,labels = mydata3$spcode,col=mycol[match(mydata3$Species,species_list)],cex=0.7)
abline(lm(log10(mydata3$ratio)~mydata3$pc1))
abline(h=1,lty=3)
dev.off()


