# Script to retrieve results of the linked selection model and combine them into figure
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022

if (!require("plotrix"))   install.packages("plotrix", dependencies = TRUE)
if (!require("ploDescToolstrix"))   install.packages("DescTools", dependencies = TRUE)
if (!require("PropCIs"))   install.packages("PropCIs", dependencies = TRUE)

data_rec <- read.table("data/recombination/RecombinationRates_AllHordeumGenes.txt",header=T)
data_hordeum <- read.table("data/recombination/HordeumGenes.txt",header=T)
  
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
    "T_boeticum",
    "T_urartu"
  )
  
species_list_abrev <- c(
    "Abi",
    "Aca",
    "Aco",
    "Alo",
    "Amu",
    "Ase",
    "Ash",
    "Asp",
    "Ata",
    "Aum",
    "Aun",
    "Tmo",
    "Tur"
  )
  
  # Colors for figures
  mycol <- rep(NA, length(species_list))
  mycol[c(which(species_list == "Ae_mutica"))] <- "#085CF8"
  mycol[c(which(species_list == "Ae_speltoides"))] <- "#1684A6"   # "royalblue"
  mycol[c(which(species_list == "Ae_sharonensis"))] <- "#3C9E49"
  mycol[c(which(species_list == "Ae_caudata"))] <- "#65AF1E"
  mycol[c(which(species_list == "Ae_longissima"))] <- "#98BB18"
  mycol[c(which(species_list == "Ae_umbellulata"))] <- "#C7C612" 
  mycol[c(which(species_list == "Ae_bicornis"))] <- "#F3CC1D"
  mycol[c(which(species_list == "Ae_comosa"))] <- "#FDAF55"
  mycol[c(which(species_list == "T_boeticum"))] <- "#FE8F7B"
  mycol[c(which(species_list == "Ae_searsii"))] <- "#FC6A9B"
  mycol[c(which(species_list == "Ae_uniaristata"))] <- "#F64497" 
  mycol[c(which(species_list == "Ae_tauschii"))] <- "#E92952"
  mycol[c(which(species_list == "T_urartu"))] <- "#D70500"

  
  # Analysis for piS, Ncat number of categories
  
  Ncat <- 20
  
  pdf(file = "figures/main/piS_piNpiS_vs_rec.pdf", height = 8, width= 12, pointsize = 18)
  
  layout(matrix(c(1,2),ncol=2))
  
  plot(NULL,xlim=c(0.01,1.8),ylim=c(0,0.05),xlab= "Recombination rate",ylab=expression(italic(pi)[S]))
  
  taby <- seq(0.049,0,-0.0026)
  
  for(SP in c(1:13)) {
    
    # Choice of the species
    SPECIES <- species_list[SP]
    
    data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
    data_focal$V1 <- gsub("_simExt","",data_focal$V1)
    data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
    data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
    data_focal <- data_focal[,c("V1","V2","V5","V7","V9","V11","V13","V15")]
    names(data_focal) <- c("G_focal","G_hordeum","t","S","N","dNdS","dN","dS")
    # Filtering out too divergent orthologs
    THRESHOLD <- 0.35

    data_focal <- data_focal[data_focal$t<THRESHOLD,]
    FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                   "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
    PATH <- paste("data/polymorphism/",SPECIES,"/",FILE,sep="")
    data_pol <- read.table(PATH,header=T)
    # To suppress the species name in some variable names 
    names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
    data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
    data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
    data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
    
    mydata <- merge(data_rec,data_hordeum,by="Gene")
    mydata <- merge(mydata,data_focal,by.x="Gene",by.y="G_hordeum")
    mydata <- merge(mydata,data_pol,by.x="G_focal",by.y="Contig_name")
    mydata <- mydata[-which(names(mydata)==c("Chromosome.y")| names(mydata)==c("Start.y"))]
    names(mydata)[which(names(mydata)=="Chromosome.x")] <- "Chromosome"
    names(mydata)[which(names(mydata)=="Start.x")] <- "Start"
    mydata$Recombination <- ifelse(mydata$Recombination<0,0,mydata$Recombination)
    mydata$CatRec  <- cut(mydata$Recombination,breaks = unique(quantile(mydata$Recombination, c(0:Ncat)/Ncat,na.rm=T)))
    mydata$PIS <- mydata$piS*mydata$nb_complete_site
    mydata$PIN <- mydata$piN*mydata$nb_complete_site
    mydata$Ldiv <- mydata$S
    mydata$SIM <- mydata$dS*mydata$Ldiv
    mydataWindow <- aggregate(mydata,by=list(mydata$CatRec),FUN = function(x) mean(x,na.rm=T))
    mydataWindow$piSyn <- mydataWindow$PIS/mydataWindow$nb_complete_site
    mydataWindow$piNonSyn <- mydataWindow$PIN/mydataWindow$nb_complete_site
    mydataWindow$f0 <- mydataWindow$piNonSyn/mydataWindow$piSyn
    mydataWindow$div <- mydataWindow$S*mydataWindow$dS/mydataWindow$S
    mydataWindow$Ne <- 12000000*mydataWindow$piSyn/mydataWindow$div # 12000000 is the estimated divergence with Hordeum
    mydataWindow$reldiv <- mydataWindow$div/mean(mydataWindow$div)
    mydataWindow$piSyn_cor <- mydataWindow$piSyn/mydataWindow$reldiv
 
    points(mydataWindow$piSyn~mydataWindow$Recombination,pch=16,col=mycol[SP])
    points(mydataWindow$piSyn~mydataWindow$Recombination)
    lines(loess(mydataWindow$piSyn~mydataWindow$Recombination),conf.level = F,col=mycol[SP])
    
    ratio <-  max(mydataWindow$piSyn,na.rm = T)/(min(mydataWindow[mydataWindow$piSyn>0,]$piSyn,na.rm = T))
    
    text(x=1.45,y=taby[SP],pos = 4, col=mycol[SP],labels = paste(species_list_abrev[SP],round(ratio,2),sep=" : "),cex=0.7)
    
  }
  
  
  plot(NULL,xlim=c(0,1.8),ylim=c(0,0.6),xlab="Recombination rate",ylab=expression(italic(pi)[N]/italic(pi)[S]))
  
  taby <- seq(0.58,0,-0.03)
  
  for(SP in c(1:13)) {
    
    # Choice of the species
    SPECIES <- species_list[SP]
    
    data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
    data_focal$V1 <- gsub("_simExt","",data_focal$V1)
    data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
    data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
    data_focal <- data_focal[,c("V1","V2","V5","V7","V9","V11","V13","V15")]
    names(data_focal) <- c("G_focal","G_hordeum","t","S","N","dNdS","dN","dS")
    # Filtering out too divergent orthologs
    THRESHOLD <- 0.35
    
    data_focal <- data_focal[data_focal$t<THRESHOLD,]
    FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                   "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
    
    PATH <- paste("data/polymorphism/",SPECIES,"/",FILE,sep="")
    data_pol <- read.table(PATH,header=T)
    # To suppress the species name in some variable names 
    names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
    data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
    data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
    data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
    
    mydata <- merge(data_rec,data_hordeum,by="Gene")
    mydata <- merge(mydata,data_focal,by.x="Gene",by.y="G_hordeum")
    mydata <- merge(mydata,data_pol,by.x="G_focal",by.y="Contig_name")
    mydata <- mydata[-which(names(mydata)==c("Chromosome.y")| names(mydata)==c("Start.y"))]
    names(mydata)[which(names(mydata)=="Chromosome.x")] <- "Chromosome"
    names(mydata)[which(names(mydata)=="Start.x")] <- "Start"
    mydata$Recombination <- ifelse(mydata$Recombination<0,0,mydata$Recombination)
    mydata$CatRec  <- cut(mydata$Recombination,breaks = unique(quantile(mydata$Recombination, c(0:Ncat)/Ncat,na.rm=T)))
    mydata$PIS <- mydata$piS*mydata$nb_complete_site
    mydata$PIN <- mydata$piN*mydata$nb_complete_site
    mydata$Ldiv <- mydata$S
    mydata$SIM <- mydata$dS*mydata$Ldiv
    mydataWindow <- aggregate(mydata,by=list(mydata$CatRec),FUN = function(x) mean(x,na.rm=T))
    mydataWindow$piSyn <- mydataWindow$PIS/mydataWindow$nb_complete_site
    mydataWindow$piNonSyn <- mydataWindow$PIN/mydataWindow$nb_complete_site
    mydataWindow$f0 <- mydataWindow$piNonSyn/mydataWindow$piSyn
    mydataWindow$div <- mydataWindow$S*mydataWindow$dS/mydataWindow$S
    mydataWindow$Ne <- 12000000*mydataWindow$piSyn/mydataWindow$div
    mydataWindow$reldiv <- mydataWindow$div/mean(mydataWindow$div)
    mydataWindow$piSyn_cor <- mydataWindow$piSyn/mydataWindow$reldiv
    
    points(mydataWindow$f0~mydataWindow$Recombination,pch=16,col=mycol[SP])
    points(mydataWindow$f0~mydataWindow$Recombination)
    lines(loess(mydataWindow$f0~mydataWindow$Recombination),conf.level = F,col=mycol[SP])
    ratio <-  max(mydataWindow$f0,na.rm = T)/(min(mydataWindow[mydataWindow$f0>0,]$f0,na.rm = T))
    text(x=1.45,y=taby[SP],pos = 4, col=mycol[SP],labels = paste(species_list_abrev[SP],round(ratio,2),sep=" : "),cex=0.7)
  }
  
  dev.off()
  
  
  
  #pdf(file = "figures/sup_mat/Ds_vs_rec.pdf", height = 8, width = 8, pointsize = 18)
  
  jpeg(file = "figures/sup_mat/Ds_vs_rec.jpeg", height = 800, width = 800, pointsize = 18)
  plot(NULL,xlim=c(0.01,1.5),ylim=c(0.7,1.4),xlab= "Recombination rate",ylab=expression("Relative "*italic(D)[S]))
  
  taby <- seq(1.38,1.05,-0.025)
  
  for(SP in c(1:13)) {
    
    # Choice of the species
    SPECIES <- species_list[SP]
    
    data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
    data_focal$V1 <- gsub("_simExt","",data_focal$V1)
    data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
    data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
    data_focal <- data_focal[,c("V1","V2","V5","V7","V9","V11","V13","V15")]
    names(data_focal) <- c("G_focal","G_hordeum","t","S","N","dNdS","dN","dS")
    # Filtering out too divergent orthologs
    THRESHOLD <- 0.35

    data_focal <- data_focal[data_focal$t<THRESHOLD,]
    FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                   "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
    
    PATH <- paste("data/polymorphism/",SPECIES,"/",FILE,sep="")
    data_pol <- read.table(PATH,header=T)
    # To suppress the species name in some variable names 
    names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
    data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
    data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
    data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
    
    mydata <- merge(data_rec,data_hordeum,by="Gene")
    mydata <- merge(mydata,data_focal,by.x="Gene",by.y="G_hordeum")
    mydata <- merge(mydata,data_pol,by.x="G_focal",by.y="Contig_name")
    mydata <- mydata[-which(names(mydata)==c("Chromosome.y")| names(mydata)==c("Start.y"))]
    names(mydata)[which(names(mydata)=="Chromosome.x")] <- "Chromosome"
    names(mydata)[which(names(mydata)=="Start.x")] <- "Start"
    mydata$Recombination <- ifelse(mydata$Recombination<0,0,mydata$Recombination)
    mydata$CatRec  <- cut(mydata$Recombination,breaks = unique(quantile(mydata$Recombination, c(0:Ncat)/Ncat,na.rm=T)))
    mydata$PIS <- mydata$piS*mydata$nb_complete_site
    mydata$PIN <- mydata$piN*mydata$nb_complete_site
    mydata$Ldiv <- mydata$S
    mydata$SIM <- mydata$dS*mydata$Ldiv
    mydataWindow <- aggregate(mydata,by=list(mydata$CatRec),FUN = function(x) mean(x,na.rm=T))
    mydataWindow$piSyn <- mydataWindow$PIS/mydataWindow$nb_complete_site
    mydataWindow$piNonSyn <- mydataWindow$PIN/mydataWindow$nb_complete_site
    mydataWindow$f0 <- mydataWindow$piNonSyn/mydataWindow$piSyn
    mydataWindow$div <- mydataWindow$S*mydataWindow$dS/mydataWindow$S
    mydataWindow$Ne <- 12000000*mydataWindow$piSyn/mydataWindow$div
    mydataWindow$reldiv <- mydataWindow$div/mean(mydataWindow$div)
    mydataWindow$piSyn_cor <- mydataWindow$piSyn/mydataWindow$reldiv
    
    points(mydataWindow$reldiv~mydataWindow$Recombination,pch=16,col=mycol[SP])
    points(mydataWindow$reldiv~mydataWindow$Recombination)
    lines(loess(mydataWindow$reldiv~mydataWindow$Recombination),conf.level = F,col=mycol[SP])
    abline(h=1,lty=2)
    ratio <-  max(mydataWindow$reldiv,na.rm = T)/(min(mydataWindow[mydataWindow$reldiv>0,]$reldiv,na.rm = T))
    text(x=0,y=taby[SP],pos = 4, col=mycol[SP],labels = paste(species_list_abrev[SP],round(ratio,2),sep=" : "),cex=0.7)
    
  }

  dev.off()
  