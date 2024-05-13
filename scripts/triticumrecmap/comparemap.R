# Comparison of recombination maps of Hordeum and Triticum
# Sylvain Glémin (CNRS Rennes, France)
# May 2024



data_recH <- read.table(file ="data/recombination/RecombinationRates_AllHordeumGenes.txt",header=T)
data_recT <- read.table(file ="data/recombination/RecombinationRates_AllTriticumGenes.txt",header=T)
data_recT$recRate <- data_recT$recRate*10^6
mareymapH <- read.table("data/recombination/MareyMap_Hordeum_SNPJHutton2012_On_EnsemblPlantsVersionHv_IBSC_PGSB_v21.txt",header=T)
mareymapT <- read.table("data/recombination/MareyMap_Triticum_aestivum_GutierrezGonzalez2019.txt",header=T)
data_hordeum <- read.table("data/recombination/HordeumGenes.txt",header=T)
data_triticum <- read.table("data/recombination/TriticumGenes.txt",header=T)

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

SPECIES <- species_list[12]

data_focalH <- read.table(paste("outputs/orthology/hordeum/",SPECIES,"_Hordeum_reciprocalbestblast.txt",sep=""),header = T)
data_focalTA <- read.table(paste("outputs/orthology/triticum/",SPECIES,"_TaestivumA_reciprocalbestblast.txt",sep=""),header = T)
data_focalTB <- read.table(paste("outputs/orthology/triticum/",SPECIES,"_TaestivumB_reciprocalbestblast.txt",sep=""),header = T)
data_focalTD <- read.table(paste("outputs/orthology/triticum/",SPECIES,"_TaestivumD_reciprocalbestblast.txt",sep=""),header = T)


df1 <- merge(data_focalH[,c(1,2)],data_focalTA[,c(1,2)],by = "Query",all = T)
df2 <- merge(df1,data_focalTB[,c(1,2)],by = "Query",all = T)
df3 <- merge(df2,data_focalTD[,c(1,2)],by = "Query",all = T)
names(df3) <- c("Gene_focal","Gene_hordeum","Gene_triticum_A","Gene_triticum_B","Gene_triticum_D")
df4 <- merge(df3,data_recH[,c(1,4)],by.x = "Gene_hordeum",by.y = "Gene",all.x=T)
names(df4)[6] <- "Rec_H"
df5 <- merge(df4,data_recT[substr(data_recT$map,2,2)=="A",c(1,6)],by.x = "Gene_triticum_A",by.y = "gene",all.x=T)
names(df5)[7] <- "Rec_A"
df6 <- merge(df5,data_recT[substr(data_recT$map,2,2)=="B",c(1,6)],by.x = "Gene_triticum_B",by.y = "gene",all.x=T)
names(df6)[8] <- "Rec_B"
data_tot <- merge(df6,data_recT[substr(data_recT$map,2,2)=="D",c(1,6)],by.x = "Gene_triticum_D",by.y = "gene",all.x=T)
names(data_tot)[9] <- "Rec_D"


hist(data_tot$Rec_H,breaks = 100)
hist(data_tot$Rec_A,breaks = 100)
hist(data_tot$Rec_B,breaks = 100)
hist(data_tot$Rec_D,breaks = 100)

mrH <- median(data_tot$Rec_H,na.rm=T)
mrA <- median(data_tot$Rec_A,na.rm=T)
mrB <- median(data_tot$Rec_B,na.rm=T)
mrD <- median(data_tot$Rec_D,na.rm=T)

(length(which(data_tot$Rec_H<=mrH & data_tot$Rec_A<=mrA))+
length(which(data_tot$Rec_H>mrH & data_tot$Rec_A>mrA)))/
  length(which(!is.na(data_tot$Rec_H) & !is.na(data_tot$Rec_A)))

(length(which(data_tot$Rec_H<=mrH & data_tot$Rec_B<=mrB))+
    length(which(data_tot$Rec_H>mrH & data_tot$Rec_B>mrB)))/
  length(which(!is.na(data_tot$Rec_H) & !is.na(data_tot$Rec_B)))

(length(which(data_tot$Rec_A<=mrA & data_tot$Rec_B<=mrB))+
    length(which(data_tot$Rec_A>mrA & data_tot$Rec_B>mrB)))/
  length(which(!is.na(data_tot$Rec_A) & !is.na(data_tot$Rec_B)))



