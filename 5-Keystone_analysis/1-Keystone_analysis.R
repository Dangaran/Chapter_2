##%######################################################%##
#                                                          #
####                     Keystones                      ####
#                                                          #
##%######################################################%##
# degree vs closeness
test <- data.frame(degree = degree(bact_network_opw_igraph)/max(degree(bact_network_opw_igraph)), closen = centr_clo(bact_network_opw_igraph)$res) 
library(ggplot2)
ggplot(test, aes(degree, closen)) + geom_point() + stat_smooth(method = "lm", col = "red")


# Keystones review
# Bacteria
# Forest
degree_illumina_for <- data.frame(ID = V(bact_network_for_igraph)$OTU_ID, degree = degree(bact_network_for_igraph)/max(degree(bact_network_for_igraph))) 
degree_illumina_for <- degree_illumina_for[order(-degree_illumina_for$degree),]
head(degree_illumina_for, 25)
degree_illumina_for[1:round(nrow(degree_illumina_opw)*0.1, 0),]

# Dehesa
degree_illumina_deh <- data.frame(ID = V(bact_network_deh_igraph)$OTU_ID, degree = degree(bact_network_deh_igraph)/max(degree(bact_network_deh_igraph))) 
degree_illumina_deh <- degree_illumina_deh[order(-degree_illumina_deh$degree),]
degree_illumina_deh

# Open woodland
degree_illumina_opw <- data.frame(ID = V(bact_network_opw_igraph)$OTU_ID, degree = degree(bact_network_opw_igraph)/max(degree(bact_network_opw_igraph))) 
degree_illumina_opw <- degree_illumina_opw[order(-degree_illumina_opw$degree),]
degree_illumina_opw


# Fungal
# Forest
degree_mycota_for <- data.frame(ID = V(fungi_network_for_igraph)$OTU_ID, degree = degree(fungi_network_for_igraph)/max(degree(fungi_network_for_igraph))) 
degree_mycota_for <- degree_mycota_for[order(-degree_mycota_for$degree),]
degree_mycota_for

# Dehesa
degree_mycota_deh <- data.frame(ID = V(fungi_network_deh_igraph)$OTU_ID, degree = degree(fungi_network_deh_igraph)/max(degree(fungi_network_deh_igraph))) 
degree_mycota_deh <- degree_mycota_deh[order(-degree_mycota_deh$degree),]
degree_mycota_deh

# Open woodland
degree_mycota_opw <- data.frame(ID = V(fungi_network_opw_igraph)$OTU_ID, degree = degree(fungi_network_opw_igraph)/max(degree(fungi_network_opw_igraph))) 
degree_mycota_opw <- degree_mycota_opw[order(-degree_mycota_opw$degree),]
degree_mycota_opw

# Degree histogram
par(mfrow=c(2,3))
hist(degree_illumina_for$degree, main = "Degree BACTERIA-Forest")
hist(degree_illumina_deh$degree, main = "Degree BACTERIA-Dehesa")
hist(degree_illumina_opw$degree, main = "Degree BACTERIA-Open woodland")

hist(degree_mycota_for$degree, main = "Degree Fungi-Forest")
hist(degree_mycota_deh$degree, main = "Degree Fungi-Dehesa")
hist(degree_mycota_opw$degree, main = "Degree Fungi-Open woodland")
dev.off()

# Keystones Degree
Bact_for_deg <- degree_illumina_for[which(degree_illumina_for$degree >= 0.8),]
Bact_deh_deg <- degree_illumina_deh[which(degree_illumina_deh$degree >= 0.8),]
Bact_opw_deg <- degree_illumina_opw[which(degree_illumina_opw$degree >= 0.8),]

Fungi_for_deg <- degree_mycota_for[which(degree_mycota_for$degree >= 0.8),]
Fungi_deh_deg <- degree_mycota_deh[which(degree_mycota_deh$degree >= 0.8),]
Fungi_opw_deg <- degree_mycota_opw[which(degree_mycota_opw$degree >= 0.8),]

Bact_for_deg$Tipo_Keyst <- "Degree"
Bact_deh_deg$Tipo_Keyst <- "Degree"
Bact_opw_deg$Tipo_Keyst <- "Degree"
Fungi_for_deg$Tipo_Keyst <- "Degree"
Fungi_deh_deg$Tipo_Keyst <- "Degree"
Fungi_opw_deg$Tipo_Keyst <- "Degree"



# BTW histogram
par(mfrow=c(2,3))
hist(centr_betw_illumina_for$centr_betw, main = "Betweness BACTERIA-Forest")
hist(centr_betw_illumina_deh$centr_betw, main = "Betweness BACTERIA-Dehesa")
hist(centr_betw_illumina_opw$centr_betw, main = "Betweness BACTERIA-Open woodland")

hist(centr_betw_mycota_for$centr_betw, main = "Betweness Fungi-Forest")
hist(centr_betw_mycota_deh$centr_betw, main = "Betweness Fungi-Dehesa")
hist(centr_betw_mycota_opw$centr_betw, main = "Betweness Fungi-Open woodland")
dev.off()

# Keystones BC

Bact_for_btw <- centr_betw_illumina_for[which(centr_betw_illumina_for$centr_betw >= 0.8),]
Bact_deh_btw <- centr_betw_illumina_deh[which(centr_betw_illumina_deh$centr_betw >= 0.5),]
Bact_opw_btw <- centr_betw_illumina_opw[which(centr_betw_illumina_opw$centr_betw >= 0.5),]

Fungi_for_btw <- centr_betw_mycota_for[which(centr_betw_mycota_for$centr_betw >= 0.5),]
Fungi_deh_btw <- centr_betw_mycota_deh[which(centr_betw_mycota_deh$centr_betw >= 0.6),]
Fungi_opw_btw <- centr_betw_mycota_opw[which(centr_betw_mycota_opw$centr_betw >= 0.6),]

Bact_for_btw$Tipo_Keyst <- "BTW"
Bact_deh_btw$Tipo_Keyst <- "BTW"
Bact_opw_btw$Tipo_Keyst <- "BTW"
Fungi_for_btw$Tipo_Keyst <- "BTW"
Fungi_deh_btw$Tipo_Keyst <- "BTW"
Fungi_opw_btw$Tipo_Keyst <- "BTW"



# Keystones merged
keystones_bact_for <- rbind(Bact_for_deg[c(1,3:7)], Bact_for_btw[c(1,3:7)])
keystones_bact_deh <- rbind(Bact_deh_deg[c(1,3:7)], Bact_deh_btw[c(1,3:7)])
keystones_bact_opw <- rbind(Bact_opw_deg[c(1,3:7)], Bact_opw_btw[c(1,3:7)])

keystones_fungi_for <- rbind(Fungi_for_deg[c(1,3:8)], Fungi_for_btw[c(1,3:8)])
keystones_fungi_deh <- rbind(Fungi_deh_deg[c(1,3:8)], Fungi_deh_btw[c(1,3:8)])
keystones_fungi_opw <- rbind(Fungi_opw_deg[c(1,3:8)], Fungi_opw_btw[c(1,3:8)])

write.csv(keystones_bact_for, "./Nuevas_keystones/keystones_bact_for.csv")
write.csv(keystones_bact_deh, "./Nuevas_keystones/keystones_bact_deh.csv")
write.csv(keystones_bact_opw, "./Nuevas_keystones/keystones_bact_opw.csv")

write.csv(keystones_fungi_for, "./Nuevas_keystones/keystones_fungi_for.csv")
write.csv(keystones_fungi_deh, "./Nuevas_keystones/keystones_fungi_deh.csv")
write.csv(keystones_fungi_opw, "./Nuevas_keystones/keystones_fungi_opw.csv")

##%######################################################%##
#                                                          #
####                Relative abundance                  ####
#                                                          #
##%######################################################%##
# Bacteria
# Forest
relab_keystones_bact_for <-  transform_sample_counts(illumina_forest_final, function(x) x / sum(x))
relab_keystones_bact_for <- subset_taxa(relab_keystones_bact_for, rownames(tax_table(relab_keystones_bact_for)) %in% keystones_bact_for[,1])
relab_keystones_bact_for

# Analysis
sample_kst_bact_for <- data.frame(otu_table(relab_keystones_bact_for))
rownames(sample_kst_bact_for) == rownames(tax_table(relab_keystones_bact_for)) # Check order of rows
sample_kst_bact_for$Phyla <- data.frame(tax_table(relab_keystones_bact_for)[,1:3])$Phylum # Create Phyla column

# Average
relab_kst_bact_for <- data.frame(rowMeans(otu_table(relab_keystones_bact_for)))
colnames(relab_kst_bact_for) <- "Relative_abundance"

rownames(relab_kst_bact_for) == rownames(tax_table(relab_keystones_bact_for)) # Check order of rows
relab_kst_bact_for$Phyla <- data.frame(tax_table(relab_keystones_bact_for)[,1:3])$Phylum # Create Phyla column
relab_kst_bact_for$Class <- data.frame(tax_table(relab_keystones_bact_for)[,1:3])$Class # Create Phyla column

phylum_kst_bact_for <- aggregate(Relative_abundance ~ Phyla, relab_kst_bact_for, sum)
phylum_kst_bact_for$Relative_abundance <- phylum_kst_bact_for$Relative_abundance*100
phylum_kst_bact_for

phylum_kst_bact_for$st_error <- aggregate(Relative_abundance ~ Phyla, relab_kst_bact_for, st_error_NA)$Relative_abundance
phylum_kst_bact_for$st_error <- phylum_kst_bact_for$st_error*100
phylum_kst_bact_for


class_kst_bact_for <- aggregate(Relative_abundance ~ Class, relab_kst_bact_for, sum)
class_kst_bact_for$Relative_abundance <- class_kst_bact_for$Relative_abundance*100
class_kst_bact_for

# Dehesa
relab_keystones_bact_deh <-  transform_sample_counts(illumina_dehesa_final, function(x) x / sum(x))
relab_keystones_bact_deh <- subset_taxa(relab_keystones_bact_deh, rownames(tax_table(relab_keystones_bact_deh)) %in% keystones_bact_deh[,1])
relab_keystones_bact_deh

# Analysis
sample_kst_bact_deh <- data.frame(otu_table(relab_keystones_bact_deh))
rownames(sample_kst_bact_deh) == rownames(tax_table(relab_keystones_bact_deh)) # Check order of rows
sample_kst_bact_deh$Phyla <- data.frame(tax_table(relab_keystones_bact_deh)[,1:3])$Phylum # Create Phyla column


# Average
relab_kst_bact_deh <- data.frame(rowMeans(otu_table(relab_keystones_bact_deh)))
colnames(relab_kst_bact_deh) <- "Relative_abundance"

rownames(relab_kst_bact_deh) == rownames(tax_table(relab_keystones_bact_deh)) # Check order of rows
relab_kst_bact_deh$Phyla <- data.frame(tax_table(relab_keystones_bact_deh)[,1:3])$Phylum # Create Phyla column
relab_kst_bact_deh$Class <- data.frame(tax_table(relab_keystones_bact_deh)[,1:3])$Class # Create Phyla column

phylum_kst_bact_deh <- aggregate(Relative_abundance ~ Phyla, relab_kst_bact_deh, sum)
phylum_kst_bact_deh$Relative_abundance <- phylum_kst_bact_deh$Relative_abundance*100
phylum_kst_bact_deh
phylum_kst_bact_deh$st_error <- aggregate(Relative_abundance ~ Phyla, relab_kst_bact_deh, st_error_NA)$Relative_abundance
phylum_kst_bact_deh$st_error <- phylum_kst_bact_deh$st_error*100
phylum_kst_bact_deh


class_kst_bact_deh <- aggregate(Relative_abundance ~ Class, relab_kst_bact_deh, sum)
class_kst_bact_deh$Relative_abundance <- class_kst_bact_deh$Relative_abundance*100
class_kst_bact_deh

# Open woodland
relab_keystones_bact_opw <-  transform_sample_counts(illumina_openwood_final, function(x) x / sum(x))
relab_keystones_bact_opw <- subset_taxa(relab_keystones_bact_opw, rownames(tax_table(relab_keystones_bact_opw)) %in% keystones_bact_opw[,1])
relab_keystones_bact_opw


# Analysis
sample_kst_bact_opw <- data.frame(otu_table(relab_keystones_bact_opw))
rownames(sample_kst_bact_opw) == rownames(tax_table(relab_keystones_bact_opw)) # Check order of rows
sample_kst_bact_opw$Phyla <- data.frame(tax_table(relab_keystones_bact_opw)[,1:3])$Phylum # Create Phyla column


# Average
relab_kst_bact_opw <- data.frame(rowMeans(otu_table(relab_keystones_bact_opw)))
colnames(relab_kst_bact_opw) <- "Relative_abundance"

rownames(relab_kst_bact_opw) == rownames(tax_table(relab_keystones_bact_opw)) # Check order of rows
relab_kst_bact_opw$Phyla <- data.frame(tax_table(relab_keystones_bact_opw)[,1:3])$Phylum # Create Phyla column
relab_kst_bact_opw$Class <- data.frame(tax_table(relab_keystones_bact_opw)[,1:3])$Class # Create Phyla column

phylum_kst_bact_opw <- aggregate(Relative_abundance ~ Phyla, relab_kst_bact_opw, sum)
phylum_kst_bact_opw$Relative_abundance <- phylum_kst_bact_opw$Relative_abundance*100
phylum_kst_bact_opw
phylum_kst_bact_opw$st_error <- aggregate(Relative_abundance ~ Phyla, relab_kst_bact_opw, st_error_NA)$Relative_abundance
phylum_kst_bact_opw$st_error <- phylum_kst_bact_opw$st_error*100
phylum_kst_bact_opw

class_kst_bact_opw <- aggregate(Relative_abundance ~ Class, relab_kst_bact_opw, sum)
class_kst_bact_opw$Relative_abundance <- class_kst_bact_opw$Relative_abundance*100
class_kst_bact_opw



# Fungal
# Forest
relab_keystones_fungi_for <-  transform_sample_counts(mycota_forest_final, function(x) x / sum(x))
relab_keystones_fungi_for <- subset_taxa(relab_keystones_fungi_for, rownames(tax_table(relab_keystones_fungi_for)) %in% keystones_fungi_for[,1])
relab_keystones_fungi_for

# Analysis
sample_kst_fungi_for <- data.frame(otu_table(relab_keystones_fungi_for))
rownames(sample_kst_fungi_for) == rownames(tax_table(relab_keystones_fungi_for)) # Check order of rows
sample_kst_fungi_for$Phyla <- data.frame(tax_table(relab_keystones_fungi_for)[,1:3])$Phylum # Create Phyla column

# Average
relab_kst_fungi_for <- data.frame(rowMeans(otu_table(relab_keystones_fungi_for)))
colnames(relab_kst_fungi_for) <- "Relative_abundance"

rownames(relab_kst_fungi_for) == rownames(tax_table(relab_keystones_fungi_for)) # Check order of rows
relab_kst_fungi_for$Phyla <- data.frame(tax_table(relab_keystones_fungi_for)[,1:3])$Phylum # Create Phyla column
relab_kst_fungi_for$Class <- data.frame(tax_table(relab_keystones_fungi_for)[,1:3])$Class # Create Phyla column

phylum_kst_fungi_for <- aggregate(Relative_abundance ~ Phyla, relab_kst_fungi_for, sum)
phylum_kst_fungi_for$Relative_abundance <- phylum_kst_fungi_for$Relative_abundance*100
phylum_kst_fungi_for
phylum_kst_fungi_for$st_error <- aggregate(Relative_abundance ~ Phyla, relab_kst_fungi_for, st_error_NA)$Relative_abundance
phylum_kst_fungi_for$st_error <- phylum_kst_fungi_for$st_error*100
phylum_kst_fungi_for

class_kst_fungi_for <- aggregate(Relative_abundance ~ Class, relab_kst_fungi_for, sum)
class_kst_fungi_for$Relative_abundance <- class_kst_fungi_for$Relative_abundance*100
class_kst_fungi_for


# Dehesa
relab_keystones_fungi_deh <-  transform_sample_counts(mycota_dehesa_final, function(x) x / sum(x))
relab_keystones_fungi_deh <- subset_taxa(relab_keystones_fungi_deh, rownames(tax_table(relab_keystones_fungi_deh)) %in% keystones_fungi_deh[,1])
relab_keystones_fungi_deh

# Analysis
sample_kst_fungi_deh <- data.frame(otu_table(relab_keystones_fungi_deh))
rownames(sample_kst_fungi_deh) == rownames(tax_table(relab_keystones_fungi_deh)) # Check order of rows
sample_kst_fungi_deh$Phyla <- data.frame(tax_table(relab_keystones_fungi_deh)[,1:3])$Phylum # Create Phyla column

# Average
relab_kst_fungi_deh <- data.frame(rowMeans(otu_table(relab_keystones_fungi_deh)))
colnames(relab_kst_fungi_deh) <- "Relative_abundance"

rownames(relab_kst_fungi_deh) == rownames(tax_table(relab_keystones_fungi_deh)) # Check order of rows
relab_kst_fungi_deh$Phyla <- data.frame(tax_table(relab_keystones_fungi_deh)[,1:3])$Phylum # Create Phyla column
relab_kst_fungi_deh$Class <- data.frame(tax_table(relab_keystones_fungi_deh)[,1:3])$Class # Create Phyla column

phylum_kst_fungi_deh <- aggregate(Relative_abundance ~ Phyla, relab_kst_fungi_deh, sum)
phylum_kst_fungi_deh$Relative_abundance <- phylum_kst_fungi_deh$Relative_abundance*100
phylum_kst_fungi_deh
phylum_kst_fungi_deh$st_error <- aggregate(Relative_abundance ~ Phyla, relab_kst_fungi_deh, st_error_NA)$Relative_abundance
phylum_kst_fungi_deh$st_error <- phylum_kst_fungi_deh$st_error*100
phylum_kst_fungi_deh

class_kst_fungi_deh <- aggregate(Relative_abundance ~ Class, relab_kst_fungi_deh, sum)
class_kst_fungi_deh$Relative_abundance <- class_kst_fungi_deh$Relative_abundance*100
class_kst_fungi_deh


# Open woodland
relab_keystones_fungi_opw <-  transform_sample_counts(mycota_openwood_final, function(x) x / sum(x))
relab_keystones_fungi_opw <- subset_taxa(relab_keystones_fungi_opw, rownames(tax_table(relab_keystones_fungi_opw)) %in% keystones_fungi_opw[,1])
relab_keystones_fungi_opw

# Analysis
sample_kst_fungi_opw <- data.frame(otu_table(relab_keystones_fungi_opw))
rownames(sample_kst_fungi_opw) == rownames(tax_table(relab_keystones_fungi_opw)) # Check order of rows
sample_kst_fungi_opw$Phyla <- data.frame(tax_table(relab_keystones_fungi_opw)[,1:3])$Phylum # Create Phyla column

# Average
relab_kst_fungi_opw <- data.frame(rowMeans(otu_table(relab_keystones_fungi_opw)))
colnames(relab_kst_fungi_opw) <- "Relative_abundance"

rownames(relab_kst_fungi_opw) == rownames(tax_table(relab_keystones_fungi_opw)) # Check order of rows
relab_kst_fungi_opw$Phyla <- data.frame(tax_table(relab_keystones_fungi_opw)[,1:3])$Phylum # Create Phyla column
relab_kst_fungi_opw$Class <- data.frame(tax_table(relab_keystones_fungi_opw)[,1:3])$Class # Create Phyla column

phylum_kst_fungi_opw <- aggregate(Relative_abundance ~ Phyla, relab_kst_fungi_opw, sum)
phylum_kst_fungi_opw$Relative_abundance <- phylum_kst_fungi_opw$Relative_abundance*100
phylum_kst_fungi_opw
phylum_kst_fungi_opw$st_error <- aggregate(Relative_abundance ~ Phyla, relab_kst_fungi_opw, st_error_NA)$Relative_abundance
phylum_kst_fungi_opw$st_error <- phylum_kst_fungi_opw$st_error*100
phylum_kst_fungi_opw

class_kst_fungi_opw <- aggregate(Relative_abundance ~ Class, relab_kst_fungi_opw, sum)
class_kst_fungi_opw$Relative_abundance <- class_kst_fungi_opw$Relative_abundance*100
class_kst_fungi_opw

##%######################################################%##
#                                                          #
####                Management analyses                 ####
#                                                          #
##%######################################################%##
phylum_kst_bact_for
phylum_kst_bact_deh
phylum_kst_bact_opw

phylum_kst_fungi_for
phylum_kst_fungi_deh
phylum_kst_fungi_opw

# Bacteria
sample_kst_bact_for
sample_kst_bact_deh
sample_kst_bact_opw


Phyla_significance_Bacteria <- function(Phylum){
  phyla_for <-subset(sample_kst_bact_for, Phyla %in% c(Phylum))
  phyla_for <- data.frame(colSums(phyla_for[-ncol(phyla_for)]))
  colnames(phyla_for) <- "Forest"
  
  
  phyla_deh <-subset(sample_kst_bact_deh, Phyla %in% c(Phylum))
  phyla_deh <- data.frame(colSums(phyla_deh[-ncol(phyla_deh)]))
  colnames(phyla_deh) <- "Dehesa"
  
  phyla_opw <-subset(sample_kst_bact_opw, Phyla %in% c(Phylum))
  phyla_opw <- data.frame(colSums(phyla_opw[-ncol(phyla_opw)]))
  colnames(phyla_opw) <- "Open_woodland"
  
  
  # Merge dataframe
  library(reshape2)
  phyla_merged <- rbind(melt(phyla_for), melt(phyla_deh), melt(phyla_opw))
  colnames(phyla_merged) <- c("Phylum", "Relative_abundance")
  anova_test <- aov(Relative_abundance ~ Phylum, phyla_merged)
  
  anova_results <- summary(anova_test)
  tukey_results <- TukeyHSD(anova_test)
  
  output <- list(data_subseted=phyla_merged,anova_results=anova_results, tukey_results=tukey_results)
  output
}

Phyla_significance_Bacteria("Verrucomicrobia")

# Fungal
sample_kst_fungi_for
sample_kst_fungi_deh
sample_kst_fungi_opw


Phyla_significance_Fungal <- function(Phylum){
  phyla_for <-subset(sample_kst_fungi_for, Phyla %in% c(Phylum))
  phyla_for <- data.frame(colSums(phyla_for[-ncol(phyla_for)]))
  colnames(phyla_for) <- "Forest"
  
  
  phyla_deh <-subset(sample_kst_fungi_deh, Phyla %in% c(Phylum))
  phyla_deh <- data.frame(colSums(phyla_deh[-ncol(phyla_deh)]))
  colnames(phyla_deh) <- "Dehesa"
  
  phyla_opw <-subset(sample_kst_fungi_opw, Phyla %in% c(Phylum))
  phyla_opw <- data.frame(colSums(phyla_opw[-ncol(phyla_opw)]))
  colnames(phyla_opw) <- "Open_woodland"
  
  
  # Merge dataframe
  library(reshape2)
  phyla_merged <- rbind(melt(phyla_for), melt(phyla_deh), melt(phyla_opw))
  colnames(phyla_merged) <- c("Phylum", "Relative_abundance")
  anova_test <- aov(Relative_abundance ~ Phylum, phyla_merged)
  
  anova_results <- summary(anova_test)
  tukey_results <- TukeyHSD(anova_test)
  
  output <- list(data_subseted=phyla_merged,anova_results=anova_results, tukey_results=tukey_results)
  output
}

Phyla_significance_Fungal("p__Mucoromycota")

