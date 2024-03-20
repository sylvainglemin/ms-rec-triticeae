# Script to generate files for linked selection analysis
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022

# Choice of the window size in cM
SizecM <- 1

# Load datafiles
data_rec <- read.table("data/recombination/RecombinationRates_AllHordeumGenes.txt",header=T)
mareymap <- read.table("data/recombination/MareyMap_Hordeum_SNPJHutton2012_On_EnsemblPlantsVersionHv_IBSC_PGSB_v21.txt",header=T)
mareymap <- mareymap[mareymap$vld==T,]
data_hordeum <- read.table("data/recombination/HordeumGenes.txt",header=T)

data_hordeum <- data_hordeum[order(data_hordeum$Chromosome,data_hordeum$Start),]
## Adding genetic distance to the hordeum dataset
data_hordeum$dist <- rep(NA,length(data_hordeum$Gene))
spanlist <- rep(0.1,7)
for(i in 1:7){
  map <- mareymap[mareymap$map==paste(i,"H",sep=""),]
  chrom <- data_hordeum[data_hordeum$Chromosome==paste("chr",i,"H",sep=""),]
  x <- map$phys
  y <- map$gen
  fit <- loess(y~x,span = spanlist[i],degree = 2,control = loess.control(surface = "direct"))
  z <- predict(fit,data.frame(x=chrom$Start))
  dz <- z[-1]-z[-length(z)]
  dz <- ifelse(dz<0,0,dz)
  z <- 0
  for(j in 1:length(dz)) z <- c(z,z[length(z)]+dz[j])
  data_hordeum[data_hordeum$Chromosome==paste("chr",i,"H",sep=""),]$dist <- z/100
}


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

fis_list <- c(0.9,
              0.6,
              0.8,
              0.7,
              0.0,
              0.9,
              0.7,
              0.0,
              0.9,
              0.8,
              0.9,
              0.85,
              0.95
)


for(i in 1:length(species_list)) {
  
# Choice of the species
SPECIES <- species_list[i]

data_focal <- read.table(paste("outputs/orthology/",SPECIES,"_Hordeum_ds.txt",sep=""),header=F,sep ="",fill=T)
data_focal$V1 <- gsub("_simExt","",data_focal$V1)
data_focal$V1 <- gsub("_lgOrf","",data_focal$V1)
data_focal$V1 <- gsub("_ESTScan","",data_focal$V1)
data_focal <- data_focal[,c("V1","V2","V5","V7","V9","V11","V13","V15")]
names(data_focal) <- c("G_focal","G_hordeum","t","S","N","dNdS","dN","dS")

# Filtering out too divergent orthologs. Inspection by eye because suggests that 0.35 is an appropriate threshold
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

data_hordeum$poscM <- SizecM*round(100*data_hordeum$dist/SizecM)/100
data_hordeumcM <- aggregate(data_hordeum[data_hordeum$Type=="protein_coding",c("Start","End")],
                            by=list("poscM"=data_hordeum[data_hordeum$Type=="protein_coding",]$poscM,"Chromosome"=data_hordeum[data_hordeum$Type=="protein_coding",]$Chromosome),
                            FUN = function(x) sum(as.numeric(x)))
data_hordeumcM$GenDens <- (data_hordeumcM$End - data_hordeumcM$Start)


## Focal species datatset 
mydata$poscM  <- SizecM*round(100*mydata$dist/SizecM)/100
mydata$PIS <- mydata$piS*mydata$nb_complete_site
mydata$PIN <- mydata$piN*mydata$nb_complete_site
mydata$Ldiv <- mydata$S
mydata$SIM <- mydata$dS*mydata$Ldiv
mydatacM <- aggregate(mydata,by=list("poscM"=mydata$poscM,"Chromosome"=mydata$Chromosome),FUN = function(x) mean(x,na.rm=T))
mydatacM <- mydatacM[-c(3:5)]
mydatacM$nb_complete_site <- aggregate(mydata$nb_complete_site,by=list("poscM"=mydata$poscM,"Chromosome"=mydata$Chromosome),FUN = function(x) sum(x,na.rm=T))$x
mydatacM$NbGene <- aggregate(mydata$Start,by=list("poscM"=mydata$poscM,"Chromosome"=mydata$Chromosome),FUN = length)$x
mydatacM$piSyn <- mydatacM$NbGene*mydatacM$PIS/mydatacM$nb_complete_site
mydatacM$piNonSyn <- mydatacM$NbGene*mydatacM$PIN/mydatacM$nb_complete_site
mydatacM$f0 <- mydatacM$piNonSyn/mydatacM$piSyn
mydatacM$div <- mydatacM$S*mydatacM$dS/mydatacM$S
mydatacM$Ne <- 12000000*mydatacM$piSyn/mydatacM$div
mydatacM$Chromosome <- gsub("H",replacement = "",x = mydatacM$Chromosome) # Suppress the H in the chromosome name


# Exporting the dataset
toExport <- mydatacM[,c("poscM","Chromosome","Start","End","Loess1","Loess2","nb_complete_site","piSyn","f0","div","Ne")]
write.table(toExport,paste("outputs/recombination/",SPECIES,SizecM,"cM.txt",sep=""),quote = F,row.names = F)
write.table(data_hordeumcM,paste("outputs/recombination/",SizecM,"cM.txt",sep=""),quote = F,row.names = F)

print(SPECIES)

}

