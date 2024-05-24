# Script to generate files for linked selection analysis
# Sylvain Glémin (CNRS Rennes, France)
# September 2021, updated November 2022

# Choice of the window size in cM
SizecM <- 1

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


for(GENOME in c("A","B","D")){  
  
  # Load datafiles
  data_rec <- read.table("data/recombination/RecombinationRates_AllTriticumGenes.txt",header=T)
  data_rec <- data_rec[substr(data_rec$map,2,2)==GENOME,]
  mareymap <- read.table("data/recombination/MareyMap_Triticum_aestivum_GutierrezGonzalez2019.txt",header=T)
  mareymap <- mareymap[mareymap$vld==T,]
  data_triticum <- read.table("data/recombination/TriticumGenes.txt",header=T)
  data_triticum <- data_triticum[substr(data_triticum$Chromosome,2,2)==GENOME,]
  
  data_triticum <- data_triticum[order(data_triticum$Chromosome,data_triticum$Start),]
  ## Adding genetic distance to the hordeum dataset
  data_triticum$dist <- rep(NA,length(data_triticum$Gene))
  spanlist <- rep(0.1,7)
  
  if(GENOME=="A") list_chrom <- c(1:4,6:7) else list_chrom <- c(1:7)
  
  for(i in list_chrom){
    map <- mareymap[mareymap$map==paste(i,GENOME,sep=""),]
    chrom <- data_triticum[data_triticum$Chromosome==paste(i,GENOME,sep=""),]
    x <- map$phys
    y <- map$gen
    fit <- loess(y~x,span = spanlist[i],degree = 2,control = loess.control(surface = "direct"))
    z <- predict(fit,data.frame(x=chrom$Start))
    dz <- z[-1]-z[-length(z)]
    dz <- ifelse(dz<0,0,dz)
    z <- 0
    for(j in 1:length(dz)) z <- c(z,z[length(z)]+dz[j])
    data_triticum[data_triticum$Chromosome==paste(i,GENOME,sep=""),]$dist <- z/100
  }
  
 
  for(i in 1:length(species_list)) {
    
    # Choice of the species
    SPECIES <- species_list[i]
    
    data_focal <- read.table(paste("outputs/orthology/triticum/",SPECIES,"_Taestivum",GENOME,"_reciprocalbestblast.txt",sep=""),header=T,sep ="",fill=T)
    data_focal$Query <- gsub("_simExt","",data_focal$Query)
    data_focal$Query <- gsub("_lgORF","",data_focal$Query)
    data_focal$Query <- gsub("_ESTScan","",data_focal$Query)
    data_focal <- data_focal[data_focal$Evalue<10^(-20),]
    data_focal <- data_focal[,c(1,2)]
    names(data_focal) <- c("Gene_focal","Gene")
    
    FILE <- ifelse(SPECIES == "Ae_tauschii" | SPECIES == "T_boeticum" | SPECIES == "T_urartu",
                   "dNdSpiNpiS_Fis_1allele","dNdSpiNpiS_Fis_output")
    PATH <- paste("data/polymorphism/",SPECIES,"/",FILE,sep="")
    data_pol <- read.table(PATH,header=T)
    # To suppress the species name in some variable names 
    names(data_pol) <- gsub(paste(SPECIES,"_",sep=""),"",names(data_pol) )
    data_pol$piS <- as.numeric(ifelse(data_pol$piS>-1,data_pol$piS,"NA"))
    data_pol$piN <- as.numeric(ifelse(data_pol$piN>-1,data_pol$piN,"NA"))
    data_pol <- data_pol[,c("Contig_name","nb_complete_site","piS","piN")]
    
    mydata <- merge(data_triticum,data_rec,by.x="Gene",by.y="gene")
    mydata <- merge(mydata,data_focal,by="Gene")
    mydata <- merge(mydata,data_pol,by.x="Gene_focal",by.y="Contig_name")
    
    data_triticum$poscM <- SizecM*round(100*data_triticum$dist/SizecM)/100
    data_triticumcM <- aggregate(data_triticum[,c("Start","End")],
                                 by=list("poscM"=data_triticum$poscM,"Chromosome"=data_triticum$Chromosome),
                                 FUN = function(x) sum(as.numeric(x)))
    data_triticumcM$GenDens <- (data_triticumcM$End - data_triticumcM$Start)
    # Adding the mean value for chrom 5A that misses a recombiantion map
    if(GENOME =="A"){
      chr5A <- data_triticum[data_triticum$Chromosome=="5A",]
      GenDens5A <- sum(chr5A$End) - sum(chr5A$Start)
      data_triticumcM <- rbind(data_triticumcM,c(NA,"5A",sum(chr5A$Start),sum(chr5A$End),GenDens5A))
    }
    
    
    ## Focal species datatset 
    mydata$poscM  <- SizecM*round(100*mydata$dist/SizecM)/100
    mydata$PIS <- mydata$piS*mydata$nb_complete_site
    mydata$PIN <- mydata$piN*mydata$nb_complete_site
    mydatacM <- aggregate(mydata,by=list("poscM"=mydata$poscM,"Chromosome"=mydata$Chromosome),FUN = function(x) mean(x,na.rm=T))
    mydatacM <- mydatacM[-c(3:5)]
    mydatacM$nb_complete_site <- aggregate(mydata$nb_complete_site,by=list("poscM"=mydata$poscM,"Chromosome"=mydata$Chromosome),FUN = function(x) sum(x,na.rm=T))$x
    mydatacM$NbGene <- aggregate(mydata$Start,by=list("poscM"=mydata$poscM,"Chromosome"=mydata$Chromosome),FUN = length)$x
    mydatacM$piSyn <- mydatacM$NbGene*mydatacM$PIS/mydatacM$nb_complete_site
    mydatacM$piNonSyn <- mydatacM$NbGene*mydatacM$PIN/mydatacM$nb_complete_site
    mydatacM$f0 <- mydatacM$piNonSyn/mydatacM$piSyn
    
    # Exporting the dataset
    toExport <- mydatacM[,c("poscM","Chromosome","Start","End","recRate","nb_complete_site","piSyn","f0")]
    write.table(toExport,paste("outputs/recombination/triticum/",SPECIES,"_genome",GENOME,"_",SizecM,"cM.txt",sep=""),quote = F,row.names = F)
    write.table(data_triticumcM,paste("outputs/recombination/triticum/Triticum_genome",GENOME,"_",SizecM,"cM.txt",sep=""),quote = F,row.names = F)
    
    print(SPECIES)
    
  }
  
}
