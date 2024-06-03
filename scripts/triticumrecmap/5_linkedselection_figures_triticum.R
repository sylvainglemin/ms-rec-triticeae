# Script to compare the results of the linked selection model with the different reference genomes
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022

if (!require("ggplot2"))   install.packages("ggplot2", dependencies = TRUE)
if (!require("plotrix"))   install.packages("plotrix", dependencies = TRUE)
if (!require("dplyr"))   install.packages("dplyr", dependencies = TRUE)

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
for(FILE in list.files(path = "outputs/recombination/triticum/",pattern = "Results_BSfit") ){
  mydata <-  rbind(mydata, read.table(paste("outputs/recombination/triticum/",FILE,sep=""),header=T))
}
names(mydata)[1] <- "Genome"
mydata <- cbind("Species"=rep(row.names(mydata[c(1:13),]),3),mydata)

dataHordeum <-  read.table("outputs/recombination/hordeum/Results_BSfit_filter300_window1cM.txt",header=T)
dataHordeum$Genome <- "H"
col <- colnames(dataHordeum)
colreorderd <- c(col[1],col[length(col)],col[c(2:(length(col)-1))])
mydata <- rbind(mydata,dataHordeum[,colreorderd])

mydata_inf <- c()
mydata_sup <- c()
mydata_mean <- c()
for(FILE in list.files(path = "outputs/recombination/hordeum/",pattern = "Bootstrap") ){
  df <- read.table(paste("outputs/recombination/hordeum/",FILE,sep=""),header=T)
  inf <- apply(df[-c(1,2)],2,function(x)  sort(x)[3])
  sup <- apply(df[-c(1,2)],2,function(x)  sort(x)[97])
  moy <- apply(df[-c(1,2)],2,mean) 
  mydata_inf <- rbind(mydata_inf,inf)
  mydata_sup <- rbind(mydata_sup,sup)
  mydata_mean <- rbind(mydata_mean,moy)
}
mydata_inf <- data.frame(mydata_inf,row.names = NULL)
mydata_inf <- cbind("Species"=species_list,mydata_inf)
mydata_sup <- data.frame(mydata_sup,row.names = NULL)
mydata_sup <- cbind("Species"=species_list,mydata_sup)
mydata_mean <- data.frame(mydata_mean,row.names = NULL)
mydata_mean <- cbind("Species"=species_list,mydata_mean)

#mydata <- mydata[-6]
# for(FILE in list.files(path = "outputs/recombination/hordeum/",pattern = "Results_BSfit_FixedFis_filter300") ){
#   mydata <-  rbind(mydata, read.table(paste("outputs/recombination/hordeum/",FILE,sep=""),header=T))
# }
# Replace T_boeticum by T_monococcum so also change species list


mydata$Species <- ifelse(mydata$Species=="T_boeticum","T_monococcum",mydata$Species)
mydata$spcode <- spcode[match(mydata$Species,species_list)]
mydata$pc1 <- pc1_mean[match(mydata$Species,species_list)]
mydata$ratio <- mydata$piS_est/mydata$piS_obs
ord <- order(aggregate(mydata$ratio,by = list(mydata$Species),median)$x)
mydata$Species <- factor(mydata$Species,levels = species_list[ord])



#owngenome <- mydata[c(12,13,18,21,27:30,32,33,35:37),]


#pdf("figures/sup_mat/linked-selection_compare_maps.pdf",width = 8,height = 6,pointsize = 18)
jpeg("figures/sup_mat/linked-selection_compare_maps.jpeg",width = 1600,height = 700,pointsize = 18)

layout(matrix(c(1,2),ncol = 2))

# Figure piMax vs PC1
x <- mydata$pc1
y <- mydata$piS_est
shape <- case_when(
  mydata$Genome=="A" ~ 1,
  mydata$Genome=="B" ~ 2,
  mydata$Genome=="D" ~ 5,
  mydata$Genome=="H" ~ 16
)
plot(x,y, 
     xlab="PC1", ylab=expression(italic(pi)[max]~(log[10]~scale)), 
     xlim= c(min(x)-0.15, max(x)+0.15), 
     ylim = c(0.1*min(y),max(y)*1.5),
     log="y",
     pch=shape, cex.lab= 1.4,
     col=mycol[match(mydata$spcode,spcode)]) 
abline(lm(log10(y[which(mydata$Genome=="A" )])~x[which(mydata$Genome=="A" )]),lty=2)
abline(lm(log10(y[which(mydata$Genome=="B" )])~x[which(mydata$Genome=="B" )]),lty=3)
abline(lm(log10(y[which(mydata$Genome=="D" )])~x[which(mydata$Genome=="D" )]),lty=4)
abline(lm(log10(y[which(mydata$Genome=="H" )])~x[which(mydata$Genome=="H" )]),lwd=2)
points(x,y, 
       xlab="PC1", ylab=expression(italic(pi)[S]~(log[10]~scale)), 
       xlim= c(min(x)-0.15, max(x)+0.15), 
       ylim = c(0.1*min(y),max(y)*1.5),
       pch=shape, cex.lab= 1.4,
       col=mycol[match(mydata$spcode,spcode)]) 

# Figure piMax vs species range
polym <- as.data.frame(t(read.delim("data/polymorphism/Summary_polymorphism_Triticeae_extended.csv", row.names=1, sep=";")))

x <- rep(polym$Species_range_grid_1.0[c(1:13)],4)
y <- mydata$piS_est
shape <- case_when(
  mydata$Genome=="A" ~ 1,
  mydata$Genome=="B" ~ 2,
  mydata$Genome=="D" ~ 5,
  mydata$Genome=="H" ~ 16
)
plot(x,y, 
     xlab="Species range (in 1000 km2)", ylab=expression(italic(pi)[max]~(log[10]~scale)), 
     xlim= c(min(x)-0.15, max(x)+0.15), 
     ylim = c(0.1*min(y),max(y)*1.5),
     log="xy",
     pch=shape, cex.lab= 1.4,
     col=mycol[match(mydata$spcode,spcode)]) 
abline(lm(log10(y[which(mydata$Genome=="A" )])~log10(x[which(mydata$Genome=="A" )])),lty=2)
abline(lm(log10(y[which(mydata$Genome=="B" )])~log10(x[which(mydata$Genome=="B" )])),lty=3)
abline(lm(log10(y[which(mydata$Genome=="D" )])~log10(x[which(mydata$Genome=="D" )])),lty=4)
abline(lm(log10(y[which(mydata$Genome=="H" )])~log10(x[which(mydata$Genome=="H" )])),lwd=2)
points(x,y, 
       xlab="PC1", ylab=expression(italic(pi)[S]~(log[10]~scale)), 
       xlim= c(min(x)-0.15, max(x)+0.15), 
       ylim = c(0.1*min(y),max(y)*1.5),
       pch=shape, cex.lab= 1.4,
       col=mycol[match(mydata$spcode,spcode)]) 

dev.off()




# pdf("figures/sup_mat/linked-selection_model.pdf",width = 8,height = 6)
# par(mar=c(10,6,4,4))
# boxplot(mydata$ratio~mydata$Species,log="y",col=mycol[match(levels(mydata$Species),species_list)],xlab=NULL,ylab=NULL,
#         main=expression(pi[max]/pi[S]),las=2)
# dev.off()
# 
# mydata2 <- mydata[mydata$spcode!="Abi" & mydata$spcode!="Ash" & mydata$spcode!="Alo" & mydata$spcode!="Ase",]
# mydatat2 <- droplevels(mydata2)
# mydata3 <- aggregate(list(mydata2[,c("pc1","ratio")]),by=list(mydata2$Species),median)
# mydata3$spcode <- aggregate(list(mydata2[,c("spcode")]),by=list(mydata2$spcode),function(x) as.character(x[1]))[,2]
# names(mydata3) <- c("Species","pc1","ratio","spcode")
# 
# pdf("figures/sup_mat/linked-selection_model_goodspecies.pdf",width = 6,height = 4)
# par(mar=c(5,5,3,4))
# plot(mydata3$ratio~mydata3$pc1,
#      ylim=c(1,max(mydata3$ratio)),
#      xlab="PC1",ylab=expression(pi[max]/pi[S]),las=2,log="y")
# points(mydata3$ratio~mydata3$pc1,pch=16,cex=0.8,col=mycol[match(mydata3$Species,species_list)],
#      xlab="PC1",ylab=expression(pi[max]/pi[S]),las=2)
# text(x=mydata3$pc1,y=mydata3$ratio*0.85,labels = mydata3$spcode,col=mycol[match(mydata3$Species,species_list)],cex=0.7)
# abline(lm(log10(mydata3$ratio)~mydata3$pc1))
# abline(h=1,lty=3)
# dev.off()


