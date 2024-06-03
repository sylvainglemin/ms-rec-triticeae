# Comparison of recombination maps of Hordeum and Triticum
# Sylvain Glémin (CNRS Rennes, France)
# May 2024

# packages to install
if (!require("ggplot2"))   install.packages("ggplot2", dependencies = TRUE)

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


# Plotting the four recombination maps

#pdf("figures/sup_mat/Suppl_Figure_Sxx_RecombinationMap_Hordeum.pdf",width = 8,height = 6)

jpeg("figures/sup_mat/RecombinationMap_Hordeum.jpeg",width = 800,height = 600)
ggplot(data=data_recH,aes(x=Start,y=Recombination)) + geom_point(cex=0.1) + theme_classic() +
  facet_wrap(~Chromosome,scales = "free_x") + 
  theme(panel.background = element_rect(fill = NA, colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5)) +
  theme(title = element_text(size = 14),axis.title = element_text(size = 14), legend.position="none",legend.title = element_text(size = 14), legend.text = element_text(size = 12),strip.text = element_text(size = 12, face = "italic")) +
  xlab("Chromosome position (in bp)") + ylab("Recombination rate (in cM/Mb)")
dev.off()

jpeg("figures/sup_mat/RecombinationMap_TriticumA.jpeg",width = 800,height = 600)
ggplot(data=data_recT[substring(data_recT$map,2,2)=="A",],aes(x=start,y=recRate)) + geom_point(cex=0.1) + theme_classic() +
  facet_wrap(~map,scales = "free_x") + 
  theme(panel.background = element_rect(fill = NA, colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5)) +
  theme(title = element_text(size = 14),axis.title = element_text(size = 14), legend.position="none",legend.title = element_text(size = 14), legend.text = element_text(size = 12),strip.text = element_text(size = 12, face = "italic")) +
  xlab("Chromosome position (in bp)") + ylab("Recombination rate (in cM/Mb)")
dev.off()

jpeg("figures/sup_mat/RecombinationMap_TriticumB.jpeg",width = 800,height = 600)
ggplot(data=data_recT[substring(data_recT$map,2,2)=="B",],aes(x=start,y=recRate)) + geom_point(cex=0.1) + theme_classic() +
  facet_wrap(~map,scales = "free_x") + 
  theme(panel.background = element_rect(fill = NA, colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5)) +
  theme(title = element_text(size = 14),axis.title = element_text(size = 14), legend.position="none",legend.title = element_text(size = 14), legend.text = element_text(size = 12),strip.text = element_text(size = 12, face = "italic")) +
  xlab("Chromosome position (in bp)") + ylab("Recombination rate (in cM/Mb)")
dev.off()

jpeg("figures/sup_mat/RecombinationMap_TriticumD.jpeg",width = 800,height = 600)
ggplot(data=data_recT[substring(data_recT$map,2,2)=="D",],aes(x=start,y=recRate)) + geom_point(cex=0.1) + theme_classic() +
  facet_wrap(~map,scales = "free_x") + 
  theme(panel.background = element_rect(fill = NA, colour = NA),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text=element_text(size = 12), axis.text.x=element_text(face= "italic",angle=90,hjust=1,vjust=0.5)) +
  theme(title = element_text(size = 14),axis.title = element_text(size = 14), legend.position="none",legend.title = element_text(size = 14), legend.text = element_text(size = 12),strip.text = element_text(size = 12, face = "italic")) +
  xlab("Chromosome position (in bp)") + ylab("Recombination rate (in cM/Mb)")
dev.off()




result <- c()

for(SPECIES in species_list) {
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
  
  data_tot$chrH <- as.numeric(substring(data_tot$Gene_hordeum,6,6))
  data_tot$chrA <- as.numeric(substring(data_tot$Gene_triticum_A,8,8))
  data_tot$chrB <- as.numeric(substring(data_tot$Gene_triticum_B,8,8))
  data_tot$chrD <- as.numeric(substring(data_tot$Gene_triticum_D,8,8))
  
  mrecH <- median(data_tot$Rec_H,na.rm=T)
  mrecA <- median(data_tot$Rec_A,na.rm=T)
  mrecB <- median(data_tot$Rec_B,na.rm=T)
  mrecD <- median(data_tot$Rec_D,na.rm=T)
  
  data_tot$catRecH <- ifelse(data_tot$Rec_H<=mrecH,"low","high")
  data_tot$catRecA <- ifelse(data_tot$Rec_A<=mrecA,"low","high")
  data_tot$catRecB <- ifelse(data_tot$Rec_B<=mrecB,"low","high")
  data_tot$catRecD <- ifelse(data_tot$Rec_D<=mrecD,"low","high")
  
  
  
  compare_pos <- data.frame(table(data_tot$chrH,data_tot$chrA,data_tot$chrB,data_tot$chrD))
  names(compare_pos) <- c("H","A","B","D","Nb")
  compare_pos <- droplevels(compare_pos[compare_pos$H!="0",])
  synteny <- compare_pos[compare_pos$H==compare_pos$A &
                           compare_pos$A==compare_pos$B &
                           compare_pos$B==compare_pos$D,]
  synteny2 <- compare_pos[compare_pos$H==compare_pos$A |
                           compare_pos$H==compare_pos$B |
                           compare_pos$H==compare_pos$D,]

  compare_rec <- data.frame(table(data_tot$catRecH,data_tot$catRecA,data_tot$catRecB,data_tot$catRecD))
  names(compare_rec) <- c("H","A","B","D","Nb")
  simrec <- compare_rec[compare_rec$H==compare_rec$A &
                          compare_rec$A==compare_rec$B &
                          compare_rec$B==compare_rec$D,]
  simrec2 <- compare_rec[compare_rec$H==compare_rec$A |
                          compare_rec$H==compare_rec$B |
                          compare_rec$H==compare_rec$D,]
  
  result <- rbind(result,
                  c(SPECIES,
                    sum(synteny$Nb)/sum(compare_pos$Nb),
                    sum(simrec$Nb)/sum(compare_rec$Nb),
                    sum(synteny2$Nb)/sum(compare_pos$Nb),
                    sum(simrec2$Nb)/sum(compare_rec$Nb)
                    )
                  )
  
}

result <- data.frame(result)
names(result) <- c("Species","Same_chr4","Same_rec4","Same_rec2","Same_rec2")

write.table(result,"outputs/recombination/compare_maps.txt",quote = F,row.names = F)

