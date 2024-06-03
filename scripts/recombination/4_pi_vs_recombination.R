# Analysis of the relationship between recombination rate and polymorphism patterns
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022


# Loading and preparing datasets ####
data_rec <- read.table(file ="data/recombination/RecombinationRates_AllHordeumGenes.txt",header=T)
mareymap <- read.table("data/recombination/MareyMap_Hordeum_SNPJHutton2012_On_EnsemblPlantsVersionHv_IBSC_PGSB_v21.txt",header=T)
# Removing invalid markers from the Marey map
mareymap <- mareymap[mareymap$vld==T,]
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



# Analysis for piS Ncat number of categories ####

Ncat <- 100
result <- c()
for(SPECIES in species_list) {
  
  # Choice of the species
  # SPECIES <- species_list[5]
  
  if(SPECIES=="T_boeticum") {
    data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
    data_focal$V1 <- gsub("_simExt","",data_focal$V1)
    data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
    data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
    data_focal <- data_focal[,c("V1","V2")]
    names(data_focal) <- c("G_focal","G_hordeum")
    data_focal$t <- 0.1
    data_focal$S <- 300
    data_focal$N <- 700
    data_focal$dNdS <- 0.1
    data_focal$dN <- 0.01
    data_focal$dS <- 0.1
  } else {
    data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
    data_focal$V1 <- gsub("_simExt","",data_focal$V1)
    data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
    data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
    data_focal <- data_focal[,c("V1","V2","V5","V7","V9","V11","V13","V15")]
    names(data_focal) <- c("G_focal","G_hordeum","t","S","N","dNdS","dN","dS")
  }
  
  # Filtering out too divergent orthologs
  # First check by eye for each dataset.
  #plot(log(data_focal$t,10)~data_focal$S)
  #abline(h=THRESHOLD,col="red")
  #hist(log(data_focal$t,10),breaks=100)
  #abline(v=THRESHOLD,col="red")
  # 0.35 works well to remove the second mode of the Ds distribution
  THRESHOLD <- 0.35
  data_focal <- data_focal[data_focal$t<THRESHOLD,]
  
  # Loading polymorphism data
  FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                 "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
  data_pol <- read.table(paste("data/polymorphism/",SPECIES,"/",FILE,sep=""),header=T)
  # To suppress the species name in some variable names 
  names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
  data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
  data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
  data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
  # Final dataset
  mydata <- merge(data_rec,data_hordeum,by="Gene")
  mydata <- merge(mydata,data_focal,by.x="Gene",by.y="G_hordeum")
  mydata <- merge(mydata,data_pol,by.x="G_focal",by.y="Contig_name")
  mydata <- mydata[-which(names(mydata)==c("Chromosome.y")| names(mydata)==c("Start.y"))]
  names(mydata)[which(names(mydata)=="Chromosome.x")] <- "Chromosome"
  names(mydata)[which(names(mydata)=="Start.x")] <- "Start"
  mydata$Recombination <- ifelse(mydata$Recombination<0,0,mydata$Recombination)
  mydata$CatRec  <- cut(mydata$Recombination,breaks = unique(quantile(mydata$Recombination, c(0:Ncat)/Ncat,na.rm=T)))
  length(levels(mydata$CatRec))
  #mydata$CatRec <- round(mydata$Loess2,2)
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

  # Fit model on piS
  filter <- which(!is.na(mydataWindow$piSyn_cor) &!is.infinite(mydataWindow$piSyn_cor) & mydataWindow$piSyn_cor>0)
  x <- mydataWindow$Recombination[filter]
  y <- mydataWindow$piSyn_cor[filter]
  n <- mydataWindow$nb_complete_site[filter]
  binit <- ifelse(mean(mydataWindow$piSyn_cor)>0.03,0.01,0.001)
  ainit <- ifelse(mean(mydataWindow$piSyn_cor)>0.03,0.01,0.001)
  w <- n/(y*(1-y))
  model <- nls(y ~ a/(1+b*exp(-(c*x))),weights = w,start = list(a=ainit,b=binit,c=1.05),algorithm = "port",lower=c(0,0,0),upper=c(100,100,100))
  coeff <- summary(model)$parameters[,1]
  sd <- summary(model)$parameters[,2]
  pimax <- coeff[1]
  pimax.sd <- sd[1]
  pimin <- coeff[1]/(1+coeff[2])
  rec.pi <- coeff[1]*coeff[2]*coeff[3]/(1+coeff[2])^2
  slope.pi <- coeff[1]*coeff[3]/4
  medianpi <- median(y)
  fit.pi <- function(x) coeff[1]/(1+coeff[2]*exp(-(coeff[3]*x)))
  pdf(paste("figures/sup_mat/",SPECIES,"_Ncat=",Ncat,"_piS.pdf",sep=""),width = 8,height = 6)
    plot(c(-100:2000)/1000,fit.pi(c(-100:2000)/1000),type="l",lwd=2,col="red",ylim=c(0,1.2*max(mydataWindow[filter,]$piS)),
      main=SPECIES,xlab="Recombination (cM/Mb)",ylab=expression(pi[S]),cex.lab=1.2)
    points(x,y,pch=16,col="grey")
    points(x,y)
  dev.off()

  # Fit model on piN/piS
  filter <- which(!is.na(mydataWindow$f0) &!is.infinite(mydataWindow$f0))
  x <- mydataWindow$Recombination[filter]
  y <- mydataWindow$f0[filter]
  n <- mydataWindow$nb_complete_site[filter]
  w <- n/(y*(1-y))
  model <- nls(y ~ a + b*exp(-c*x),start = list(a=0.1,b=0.1,c=0.1),algorithm = "port",lower=c(-1,0,0),upper=c(2,2,2))
  coeff <- summary(model)$parameters[,1]
  sd <- summary(model)$parameters[,2]
  f0min <- coeff[1]
  f0min.sd <- sd[1]
  f0max <- coeff[1]+coeff[2]
  rec.f0 <- coeff[3]
  medianf0 <- median(y)
  fit.f0 <- function(x) coeff[1]+coeff[2]*exp(-(coeff[3]*x))
  pdf(paste("figures/sup_mat/",SPECIES,"_Ncat=",Ncat,"_piNpiS.pdf",sep=""),width = 8,height = 6)
    plot(c(0:2000)/1000,fit.f0(c(0:2000)/1000),type="l",lwd=2,col="red",ylim=c(0,1.1*max(mydataWindow[filter,]$f0)),
      main=SPECIES,xlab="Recombination (cM/Mb)",ylab=expression(pi[N]/pi[S]))
    points(x,y,pch=16,col="grey")
    points(x,y)
  dev.off()
  
  stat <- data.frame(t(c(pimin,pimax,rec.pi,slope.pi,medianpi,f0min,f0max,rec.f0,medianf0,pimax.sd,f0min.sd)))
  names(stat) <- c("piSmin","piSmax","rec.pi","slope.pi","medianpiS","f0min","f0max","rec.f0","medianf0","pimax.sd","f0min.sd")
  stat$species <- SPECIES
  result <-rbind(result,stat)
  write.table(result,paste("outputs/recombination/hordeum/stat_recombination_per_species_Ncat=",Ncat,".txt",sep=""),sep = "\t",quote = F,row.names = F)
  print(SPECIES)
}



# Analysis for piN/piS: 2 categories ####

Ncat <- 2

result <- c()

for(SP in species_list) {
  
  # Choice of the species
  #SPECIES <- species_list[8]
  SPECIES <- SP
  
  if(SPECIES=="T_boeticum") {
    data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
    data_focal$V1 <- gsub("_simExt","",data_focal$V1)
    data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
    data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
    data_focal <- data_focal[,c("V1","V2")]
    names(data_focal) <- c("G_focal","G_hordeum")
    data_focal$t <- 0.1
    data_focal$S <- 300
    data_focal$N <- 700
    data_focal$dNdS <- 0.1
    data_focal$dN <- 0.01
    data_focal$dS <- 0.1
  } else {
    data_focal <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
    data_focal$V1 <- gsub("_simExt","",data_focal$V1)
    data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
    data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
    data_focal <- data_focal[,c("V1","V2","V5","V7","V9","V11","V13","V15")]
    names(data_focal) <- c("G_focal","G_hordeum","t","S","N","dNdS","dN","dS")
  }
  
  # Filtering out too divergent orthologs
  # First check by eye for each dataset.
  #plot(log(data_focal$t,10)~data_focal$S)
  #abline(h=THRESHOLD,col="red")
  #hist(log(data_focal$t,10),breaks=100)
  #abline(v=THRESHOLD,col="red")
  # 0.35 works well to remove the second mode of the Ds distribution
  THRESHOLD <- 0.35
  data_focal <- data_focal[data_focal$t<THRESHOLD,]
  
  # Loading polymorphism data
  FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                 "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
  data_pol <- read.table(paste("data/polymorphism/",SPECIES,"/",FILE,sep=""),header=T)
  # To suppress the species name in some variable names 
  names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
  data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
  data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
  data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
  # Final dataset
  mydata <- merge(data_rec,data_hordeum,by="Gene")
  mydata <- merge(mydata,data_focal,by.x="Gene",by.y="G_hordeum")
  mydata <- merge(mydata,data_pol,by.x="G_focal",by.y="Contig_name")
  mydata <- mydata[-which(names(mydata)==c("Chromosome.y")| names(mydata)==c("Start.y"))]
  names(mydata)[which(names(mydata)=="Chromosome.x")] <- "Chromosome"
  names(mydata)[which(names(mydata)=="Start.x")] <- "Start"
  mydata$Recombination <- ifelse(mydata$Recombination<0,0,mydata$Recombination)
  mydata$CatRec  <- cut(mydata$Recombination,breaks = unique(quantile(mydata$Recombination, c(0:Ncat)/Ncat,na.rm=T)))
  length(levels(mydata$CatRec))
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
  
  stat <- data.frame(t(c(mydataWindow$Recombination,mydataWindow$f0)))
  names(stat) <- c("reclow","rechigh","f0low","f0high")
  stat$species <- SPECIES
  result <-rbind(result,stat)
  write.table(result,paste("outputs/recombination/hordeum/stat_recombination_per_species_Ncat=100.txtstat_recombination_per_species_Ncat=",Ncat,".txt",sep=""),sep = "\t",quote = F,row.names = F)
  
  print(SPECIES)
}
