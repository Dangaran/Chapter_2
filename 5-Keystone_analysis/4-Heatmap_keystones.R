keystones_bact_for_mod <- read.csv("Nuevas_keystones/keystones_bact_for_mod_mod.csv", h = T, sep = ";")
keystones_bact_deh_mod <- read.csv("Nuevas_keystones/keystones_bact_deh_mod.csv", h = T, sep = ";")
keystones_bact_opw_mod <- read.csv("Nuevas_keystones/keystones_bact_opw_mod.csv", h = T, sep = ";")

keystones_fungi_for_mod <- read.csv("Nuevas_keystones/keystones_fungi_for_mod.csv", h = T, sep = ";")
keystones_fungi_deh_mod <- read.csv("Nuevas_keystones/keystones_fungi_deh_mod.csv", h = T, sep = ";")
keystones_bact_for_mod <- read.csv("Nuevas_keystones/keystones_bact_for_mod.csv", h = T, sep = ";")


##%######################################################%##
#                                                          #
####             Heatmap Keystones_bact_forest          ####
#                                                          #
##%######################################################%##

relab_keystones_bact_for <- transform_sample_counts(illumina_forest_final, function(x) x / sum(x))
relab_keystones_bact_for <- subset_taxa(relab_keystones_bact_for, rownames(tax_table(relab_keystones_bact_for)) %in% keystones_bact_for_mod[,1])
relab_keystones_bact_for


#relab_keystones_bact_for = tax_glom(relab_keystones_bact_for, "Phylum")


keystones_bact_for_data <- data.frame(sample_data(relab_keystones_bact_for))
keystones_bact_for_data$Tree_inf <- keystones_bact_for_data$Porcentaje_defol
keystones_bact_for_data$Tree_inf <- ifelse(is.na(keystones_bact_for_data$Tree_inf) == TRUE, 150, keystones_bact_for_data$Tree_inf)
colnames(keystones_bact_for_data)[c(15,30,31, 36, 43,46,45,50,49,53)] <- c("Defol", "Temp", "Prec","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH")
sample_data(relab_keystones_bact_for) <- keystones_bact_for_data


# MicrobiomeSEQ 

taxa_names(relab_keystones_bact_for) <- paste("OTU", seq(1,length(as.vector(tax_table(relab_keystones_bact_for)[,2]))), sep = "_")
colnames(keystones_bact_for_data)
correlacion_nut_func <- taxa.env.correlation(relab_keystones_bact_for, grouping_column = "Manejo_Conteo", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                             select.variables = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"))

correlacion_nut_func <- as.data.frame(correlacion_nut_func)
# Add  taxonomy to df
relab_keystones_bact_for_tax_table <- tax_table(relab_keystones_bact_for)[,1:5]

correlacion_nut_func$Phylum <- as.factor(relab_keystones_bact_for_tax_table[match(correlacion_nut_func$Taxa, rownames(relab_keystones_bact_for_tax_table)), 2])

correlacion_nut_func$OTU_ID <- as.factor(correlacion_nut_func$Taxa)

# Change otus name to alphabetical order as in tables
correlacion_nut_func$OTU_ID <- ifelse(correlacion_nut_func$OTU_ID == "OTU_14","OTU_01", 
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_15","OTU_02",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_18","OTU_03",     
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_2","OTU_04",  
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_3","OTU_05",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_6","OTU_06",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_10","OTU_07",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_11","OTU_08",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_9","OTU_09",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_8","OTU_10",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_12","OTU_11",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_1","OTU_12",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_16","OTU_13",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_17","OTU_14",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_19","OTU_15",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_20","OTU_16",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_21","OTU_17",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_4","OTU_18",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_5","OTU_19",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_7","OTU_20",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_13","OTU_21",correlacion_nut_func$OTU_ID)))))))))))))))))))))


correlacion_nut_func$Phylum_mod <- as.factor(paste(correlacion_nut_func$Phylum, correlacion_nut_func$OTU_ID))
correlacion_nut_func$Type <- "for"

# Sort variables
correlacion_nut_func$Env <- factor(correlacion_nut_func$Env, levels=c("Temp", "Prec", "pH", "Defol", "Tree_inf", "Org_C", "NH4_N", "NO3_N", "Av_P", "RH", "Ramm", "Rnit"))

# Sort by phyla
library(forcats)
correlacion_nut_func$Phylum_mod <- fct_rev(correlacion_nut_func$Phylum_mod)



library(ggplot2)
ggplot(aes(x=Type, y=Phylum_mod, fill=Correlation), data=correlacion_nut_func)+
  geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", limits=c(-1, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)+
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),
                                                                 text=element_text(family="Times New Roman", face="plain", size=19))


##%######################################################%##
#                                                          #
####        Heatmap Keystones_bact_dehesa               ####
#                                                          #
##%######################################################%##

relab_keystones_bact_deh <- transform_sample_counts(illumina_dehesa_final, function(x) x / sum(x))
relab_keystones_bact_deh <- subset_taxa(relab_keystones_bact_deh, rownames(tax_table(relab_keystones_bact_deh)) %in% keystones_bact_deh_mod[,1])
relab_keystones_bact_deh


#relab_keystones_bact_deh = tax_glom(relab_keystones_bact_deh, "Phylum")


keystones_bact_deh_data <- data.frame(sample_data(relab_keystones_bact_deh))
keystones_bact_deh_data$Tree_inf <- keystones_bact_deh_data$Porcentaje_defol
keystones_bact_deh_data$Tree_inf <- ifelse(is.na(keystones_bact_deh_data$Tree_inf) == TRUE, 150, keystones_bact_deh_data$Tree_inf)
colnames(keystones_bact_deh_data)[c(15,30,31, 36, 43,46,45,50,49,53)] <- c("Defol", "Temp", "Prec","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH")
sample_data(relab_keystones_bact_deh) <- keystones_bact_deh_data


# MicrobiomeSEQ 

taxa_names(relab_keystones_bact_deh) <- paste("OTU", seq(1,length(as.vector(tax_table(relab_keystones_bact_deh)[,2]))), sep = "_")
colnames(keystones_bact_deh_data)
correlacion_nut_func <- taxa.env.correlation(relab_keystones_bact_deh, grouping_column = "Manejo_Conteo", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                             select.variables = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"))

correlacion_nut_func <- as.data.frame(correlacion_nut_func)
# Add  taxonomy to df
relab_keystones_bact_deh_tax_table <- tax_table(relab_keystones_bact_deh)[,1:5]

correlacion_nut_func$Phylum <- as.factor(relab_keystones_bact_deh_tax_table[match(correlacion_nut_func$Taxa, rownames(relab_keystones_bact_deh_tax_table)), 2])

correlacion_nut_func$OTU_ID <- as.factor(correlacion_nut_func$Taxa)

# Change otus name to alphabetical order as in tables
correlacion_nut_func$OTU_ID <- ifelse(correlacion_nut_func$OTU_ID == "OTU_11","OTU_01", 
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_2","OTU_02",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_3","OTU_03",     
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_4","OTU_04",  
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_5","OTU_05",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_1","OTU_06",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_10","OTU_07",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_9","OTU_08",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_12","OTU_09",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_13","OTU_10",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_6","OTU_11",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_7","OTU_12",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_8","OTU_13", correlacion_nut_func$OTU_ID)))))))))))))


correlacion_nut_func$Phylum_mod <- as.factor(paste(correlacion_nut_func$Phylum, correlacion_nut_func$OTU_ID))

# Sort variables
correlacion_nut_func$Env <- factor(correlacion_nut_func$Env, levels=c("Temp", "Prec", "pH", "Defol", "Tree_inf", "Org_C", "NH4_N", "NO3_N", "Av_P", "RH", "Ramm", "Rnit"))

# Sort by phyla
library(forcats)
correlacion_nut_func$Phylum_mod <- fct_rev(correlacion_nut_func$Phylum_mod)


correlacion_nut_func$Type <- "deh"

library(ggplot2)
ggplot(aes(x=Type, y=Phylum_mod, fill=Correlation), data=correlacion_nut_func)+
  geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", limits=c(-1, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)+
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),
                                                                 text=element_text(family="Times New Roman", face="plain", size=19))




##%######################################################%##
#                                                          #
####             Heatmap Keystones_bact_openwood        ####
#                                                          #
##%######################################################%##

relab_keystones_bact_opw <- transform_sample_counts(illumina_openwood_final, function(x) x / sum(x))
relab_keystones_bact_opw <- subset_taxa(relab_keystones_bact_opw, rownames(tax_table(relab_keystones_bact_opw)) %in% keystones_bact_opw_mod[,1])
relab_keystones_bact_opw


keystones_bact_opw_data <- data.frame(sample_data(relab_keystones_bact_opw))
keystones_bact_opw_data$Tree_inf <- keystones_bact_opw_data$Porcentaje_defol
keystones_bact_opw_data$Tree_inf <- ifelse(is.na(keystones_bact_opw_data$Tree_inf) == TRUE, 150, keystones_bact_opw_data$Tree_inf)
colnames(keystones_bact_opw_data)[c(15,30,31, 36, 43,46,45,50,49,53)] <- c("Defol", "Temp", "Prec","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH")
sample_data(relab_keystones_bact_opw) <- keystones_bact_opw_data


# MicrobiomeSEQ 

taxa_names(relab_keystones_bact_opw) <- paste("OTU", seq(1,length(as.vector(tax_table(relab_keystones_bact_opw)[,2]))), sep = "_")
colnames(keystones_bact_opw_data)
correlacion_nut_func <- taxa.env.correlation(relab_keystones_bact_opw, grouping_column = "Manejo_Conteo", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                             select.variables = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"))

correlacion_nut_func <- as.data.frame(correlacion_nut_func)
# Add  taxonomy to df
relab_keystones_bact_opw_tax_table <- tax_table(relab_keystones_bact_opw)[,1:5]

correlacion_nut_func$Phylum <- as.factor(relab_keystones_bact_opw_tax_table[match(correlacion_nut_func$Taxa, rownames(relab_keystones_bact_opw_tax_table)), 2])

correlacion_nut_func$OTU_ID <- as.factor(correlacion_nut_func$Taxa)

# Change otus name to alphabetical order as in tables
correlacion_nut_func$OTU_ID <- ifelse(correlacion_nut_func$OTU_ID == "OTU_14","OTU_01", 
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_15","OTU_02",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_8","OTU_03",     
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_16","OTU_04",  
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_5","OTU_05",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_6","OTU_06",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_10","OTU_07",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_9","OTU_08",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_11","OTU_09",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_1","OTU_10",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_17","OTU_11",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_18","OTU_12",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_19","OTU_13",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_2","OTU_14",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_3","OTU_15",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_4","OTU_16",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_7","OTU_17",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_12","OTU_18",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_13","OTU_19",correlacion_nut_func$OTU_ID)))))))))))))))))))


correlacion_nut_func$Phylum_mod <- as.factor(paste(correlacion_nut_func$Phylum, correlacion_nut_func$OTU_ID))

# Sort variables
correlacion_nut_func$Env <- factor(correlacion_nut_func$Env, levels=c("Temp", "Prec", "pH", "Defol", "Tree_inf", "Org_C", "NH4_N", "NO3_N", "Av_P", "RH", "Ramm", "Rnit"))

# Sort by phyla
library(forcats)
correlacion_nut_func$Phylum_mod <- fct_rev(correlacion_nut_func$Phylum_mod)

correlacion_nut_func$Type <- "opw"

library(ggplot2)
ggplot(aes(x=Type, y=Phylum_mod, fill=Correlation), data=correlacion_nut_func)+
  geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", limits=c(-1, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)+
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),
                                                                 text=element_text(family="Times New Roman", face="plain", size=19))



##%######################################################%##
#                                                          #
####             Heatmap Keystones_fungi_forest         ####
#                                                          #
##%######################################################%##

relab_keystones_fungi_for <- transform_sample_counts(mycota_forest_final, function(x) x / sum(x))
relab_keystones_fungi_for <- subset_taxa(relab_keystones_fungi_for, rownames(tax_table(relab_keystones_fungi_for)) %in% keystones_fungi_for_mod[,1])
relab_keystones_fungi_for


keystones_fungi_for_data <- data.frame(sample_data(relab_keystones_fungi_for))
keystones_fungi_for_data$Tree_inf <- keystones_fungi_for_data$Porcentaje_defol
keystones_fungi_for_data$Tree_inf <- ifelse(is.na(keystones_fungi_for_data$Tree_inf) == TRUE, 150, keystones_fungi_for_data$Tree_inf)
colnames(keystones_fungi_for_data)[c(15,30,31, 36, 43,46,45,50,49,53)] <- c("Defol", "Temp", "Prec","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH")
sample_data(relab_keystones_fungi_for) <- keystones_fungi_for_data


# MicrobiomeSEQ 

taxa_names(relab_keystones_fungi_for) <- paste("OTU", seq(1,length(as.vector(tax_table(relab_keystones_fungi_for)[,2]))), sep = "_")
colnames(keystones_fungi_for_data)
correlacion_nut_func <- taxa.env.correlation(relab_keystones_fungi_for, grouping_column = "Manejo_Conteo", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                             select.variables = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"))

correlacion_nut_func <- as.data.frame(correlacion_nut_func)
# Add  taxonomy to df
relab_keystones_fungi_for_tax_table <- tax_table(relab_keystones_fungi_for)[,1:5]

correlacion_nut_func$Phylum <- as.factor(relab_keystones_fungi_for_tax_table[match(correlacion_nut_func$Taxa, rownames(relab_keystones_fungi_for_tax_table)), 2])

correlacion_nut_func$OTU_ID <- as.factor(correlacion_nut_func$Taxa)

# Change otus name to alphabetical order as in tables
correlacion_nut_func$OTU_ID <- ifelse(correlacion_nut_func$OTU_ID == "OTU_10","OTU_01", 
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_11","OTU_02",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_12","OTU_03",     
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_13","OTU_04",  
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_15","OTU_05",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_16","OTU_06",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_17","OTU_07",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_18","OTU_08",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_19","OTU_09",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_2","OTU_10",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_3","OTU_11",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_4","OTU_12",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_6","OTU_13",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_7","OTU_14",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_8","OTU_15",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_9","OTU_16",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_1","OTU_17",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_14","OTU_18",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_5","OTU_19",correlacion_nut_func$OTU_ID)))))))))))))))))))


correlacion_nut_func$Phylum_mod <- as.factor(paste(correlacion_nut_func$Phylum, correlacion_nut_func$OTU_ID))

# Sort variables
correlacion_nut_func$Env <- factor(correlacion_nut_func$Env, levels=c("Temp", "Prec", "pH", "Defol", "Tree_inf", "Org_C", "NH4_N", "NO3_N", "Av_P", "RH", "Ramm", "Rnit"))

# Sort by phyla
library(forcats)
correlacion_nut_func$Phylum_mod <- fct_rev(correlacion_nut_func$Phylum_mod)

correlacion_nut_func$Type <- "for"

library(ggplot2)
ggplot(aes(x=Type, y=Phylum_mod, fill=Correlation), data=correlacion_nut_func)+
  geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", limits=c(-1, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)+
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),
                                                                 text=element_text(family="Times New Roman", face="plain", size=19))



##%######################################################%##
#                                                          #
####             Heatmap Keystones_fungi_dehesa         ####
#                                                          #
##%######################################################%##

relab_keystones_fungi_deh <- transform_sample_counts(mycota_dehesa_final, function(x) x / sum(x))
relab_keystones_fungi_deh <- subset_taxa(relab_keystones_fungi_deh, rownames(tax_table(relab_keystones_fungi_deh)) %in% keystones_fungi_deh_mod[,1])
relab_keystones_fungi_deh


keystones_fungi_deh_data <- data.frame(sample_data(relab_keystones_fungi_deh))
keystones_fungi_deh_data$Tree_inf <- keystones_fungi_deh_data$Porcentaje_defol
keystones_fungi_deh_data$Tree_inf <- ifelse(is.na(keystones_fungi_deh_data$Tree_inf) == TRUE, 150, keystones_fungi_deh_data$Tree_inf)
colnames(keystones_fungi_deh_data)[c(15,30,31, 36, 43,46,45,50,49,53)] <- c("Defol", "Temp", "Prec","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH")
sample_data(relab_keystones_fungi_deh) <- keystones_fungi_deh_data


# MicrobiomeSEQ 

taxa_names(relab_keystones_fungi_deh) <- paste("OTU", seq(1,length(as.vector(tax_table(relab_keystones_fungi_deh)[,2]))), sep = "_")
colnames(keystones_fungi_deh_data)
correlacion_nut_func <- taxa.env.correlation(relab_keystones_fungi_deh, grouping_column = "Manejo_Conteo", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                             select.variables = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"))

correlacion_nut_func <- as.data.frame(correlacion_nut_func)
# Add  taxonomy to df
relab_keystones_fungi_deh_tax_table <- tax_table(relab_keystones_fungi_deh)[,1:5]

correlacion_nut_func$Phylum <- as.factor(relab_keystones_fungi_deh_tax_table[match(correlacion_nut_func$Taxa, rownames(relab_keystones_fungi_deh_tax_table)), 2])

correlacion_nut_func$OTU_ID <- as.factor(correlacion_nut_func$Taxa)

# Change otus name to alphabetical order as in tables
correlacion_nut_func$OTU_ID <- ifelse(correlacion_nut_func$OTU_ID == "OTU_1","OTU_01", 
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_10","OTU_02",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_11","OTU_03",     
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_12","OTU_04",  
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_13","OTU_05",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_14","OTU_06",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_15","OTU_07",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_16","OTU_08",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_17","OTU_09",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_2","OTU_10",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_3","OTU_11",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_4","OTU_12",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_5","OTU_13",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_8","OTU_14",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_9","OTU_15",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_6","OTU_16",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_7","OTU_17", correlacion_nut_func$OTU_ID)))))))))))))))))


correlacion_nut_func$Phylum_mod <- as.factor(paste(correlacion_nut_func$Phylum, correlacion_nut_func$OTU_ID))

# Sort variables
correlacion_nut_func$Env <- factor(correlacion_nut_func$Env, levels=c("Temp", "Prec", "pH", "Defol", "Tree_inf", "Org_C", "NH4_N", "NO3_N", "Av_P", "RH", "Ramm", "Rnit"))

# Sort by phyla
library(forcats)
correlacion_nut_func$Phylum_mod <- fct_rev(correlacion_nut_func$Phylum_mod)

correlacion_nut_func$Type <- "deh"

library(ggplot2)
ggplot(aes(x=Type, y=Phylum_mod, fill=Correlation), data=correlacion_nut_func)+
  geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", limits=c(-1, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)+
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),
                                                                 text=element_text(family="Times New Roman", face="plain", size=19))


##%######################################################%##
#                                                          #
####             Heatmap Keystones_fungi_open woodland  ####
#                                                          #
##%######################################################%##

relab_keystones_fungi_opw <- transform_sample_counts(mycota_openwood_final, function(x) x / sum(x))
relab_keystones_fungi_opw <- subset_taxa(relab_keystones_fungi_opw, rownames(tax_table(relab_keystones_fungi_opw)) %in% keystones_fungi_opw_mod[,1])
relab_keystones_fungi_opw


keystones_fungi_opw_data <- data.frame(sample_data(relab_keystones_fungi_opw))
keystones_fungi_opw_data$Tree_inf <- keystones_fungi_opw_data$Porcentaje_defol
keystones_fungi_opw_data$Tree_inf <- ifelse(is.na(keystones_fungi_opw_data$Tree_inf) == TRUE, 150, keystones_fungi_opw_data$Tree_inf)
colnames(keystones_fungi_opw_data)[c(15,30,31, 36, 43,46,45,50,49,53)] <- c("Defol", "Temp", "Prec","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH")
sample_data(relab_keystones_fungi_opw) <- keystones_fungi_opw_data


# MicrobiomeSEQ 

taxa_names(relab_keystones_fungi_opw) <- paste("OTU", seq(1,length(as.vector(tax_table(relab_keystones_fungi_opw)[,2]))), sep = "_")
colnames(keystones_fungi_opw_data)
correlacion_nut_func <- taxa.env.correlation(relab_keystones_fungi_opw, grouping_column = "Manejo_Conteo", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                             select.variables = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"))

correlacion_nut_func <- as.data.frame(correlacion_nut_func)
# Add  taxonomy to df
relab_keystones_fungi_opw_tax_table <- tax_table(relab_keystones_fungi_opw)[,1:5]

correlacion_nut_func$Phylum <- as.factor(relab_keystones_fungi_opw_tax_table[match(correlacion_nut_func$Taxa, rownames(relab_keystones_fungi_opw_tax_table)), 2])

correlacion_nut_func$OTU_ID <- as.factor(correlacion_nut_func$Taxa)

# Change otus name to alphabetical order as in tables
correlacion_nut_func$OTU_ID <- ifelse(correlacion_nut_func$OTU_ID == "OTU_1","OTU_01", 
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_11","OTU_02",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_12","OTU_03",     
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_13","OTU_04",  
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_14","OTU_05",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_16","OTU_06",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_18","OTU_07",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_19","OTU_08",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_2","OTU_09",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_20","OTU_10",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_21","OTU_11",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_22","OTU_12",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_23","OTU_13",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_24","OTU_14",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_4","OTU_15",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_5","OTU_16",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_6","OTU_17",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_8","OTU_18",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_10","OTU_19",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_17","OTU_20",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_3","OTU_21",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_9","OTU_22",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_15","OTU_23",
                                      ifelse(correlacion_nut_func$OTU_ID == "OTU_7","OTU_24",correlacion_nut_func$OTU_ID))))))))))))))))))))))))


correlacion_nut_func$Phylum_mod <- as.factor(paste(correlacion_nut_func$Phylum, correlacion_nut_func$OTU_ID))

# Sort variables
correlacion_nut_func$Env <- factor(correlacion_nut_func$Env, levels=c("Temp", "Prec", "pH", "Defol", "Tree_inf", "Org_C", "NH4_N", "NO3_N", "Av_P", "RH", "Ramm", "Rnit"))

# Sort by phyla
library(forcats)
correlacion_nut_func$Phylum_mod <- fct_rev(correlacion_nut_func$Phylum_mod)

correlacion_nut_func$Type <- "opw"

library(ggplot2)
ggplot(aes(x=Type, y=Phylum_mod, fill=Correlation), data=correlacion_nut_func)+
  geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", limits=c(-1, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)+
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),
                                                                 text=element_text(family="Times New Roman", face="plain", size=19))

