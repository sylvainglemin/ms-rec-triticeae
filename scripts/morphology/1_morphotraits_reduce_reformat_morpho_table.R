#Part 1. Preliminar subsetting and reformatting of the morphological data"
# author: Concetta Burgarella (Uppsala University)
# date: November 2022 (Change name of T boeticum to T monococcum)

if (!require("psych"))   install.packages("psych", dependencies = TRUE)
if (!require("missMDA"))   install.packages("missMDA", dependencies = TRUE)
if (!require("FactoMineR"))   install.packages("FactoMineR", dependencies = TRUE)
if (!require("dplyr"))   install.packages("dplyr", dependencies = TRUE)


#This table include the measures directly taken on the plants (these will be imputed) and summary stats (mean, SD, other indexes) that will not be imputed and will be calculated again on the imputed data.
load("data/non_imputed_morpho_data.Rdata")

# Select species to match polymorphism estimates: include all Aegilops/Triticum, Taeniaterium caputmedusae, Hordeum spontaneum, Secale strictum
non_imputed_morpho <- rbind(subset(non_imputed_data, grepl("aeg", non_imputed_data$species)),
                            subset(non_imputed_data, grepl("tri_", non_imputed_data$species)),
                            subset(non_imputed_data, grepl("tae_", non_imputed_data$species)),
                            subset(non_imputed_data, grepl("_spo", non_imputed_data$species)),
                            subset(non_imputed_data, grepl("_str", non_imputed_data$species)))

# Recode species name for the article
non_imputed_morpho$species <- recode(non_imputed_morpho$species, aeg_bic="Abi", aeg_cau="Aca", aeg_com="Aco", aeg_lon="Alo", aeg_mut="Amu", aeg_sea="Ase", aeg_sha="Ash", aeg_spe="Asp", aeg_tau="Ata", aeg_umb="Aum", aeg_una="Aun", hor_spo="Hsp", sec_str="Str", tae_cap="Tca", tri_mon="Tmo", tri_ura="Tur")
non_imputed_morpho$ind <- gsub("aeg_bic", "Abi", non_imputed_morpho$ind) %>%
  gsub("aeg_cau", "Aca", .) %>%
  gsub("aeg_com", "Aco", .) %>%
  gsub("aeg_lon", "Alo", .) %>%
  gsub("aeg_mut", "Amu", .) %>%
  gsub("aeg_sea", "Ase", .) %>%
  gsub("aeg_sha", "Ash", .) %>%
  gsub("aeg_spe", "Asp", .) %>%
  gsub("aeg_tau", "Ata", .) %>%
  gsub("aeg_umb", "Aum", .) %>%
  gsub("aeg_una", "Aun", .) %>%
  gsub("hor_spo", "Hsp", .) %>%
  gsub("sec_str", "Str", .) %>%
  gsub("tae_cap", "Tca", .) %>%
  gsub("tri_mon", "Tmo", .) %>%
  gsub("tri_ura", "Tur", .)

# List the variables 
colnames(non_imputed_morpho)

# Choose exclusively the variables we want to include in the article (columns)
non_imputed_morpho_final <- non_imputed_morpho[,-(32:52)]

# change some columns names
names(non_imputed_morpho_final)[2] <- "species_code"
names(non_imputed_morpho_final)[3] <- "ind_code"

# Write file with morpho measures for Supplementary Table
write.csv(non_imputed_morpho_final, "outputs/sup_mat_publication/Supplementary_table_Sx_Morpho_traits_original.csv",  quote=F, row.names = F)

