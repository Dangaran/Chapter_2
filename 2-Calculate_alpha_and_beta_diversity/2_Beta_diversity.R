source('./accessory_functions/veganotu.R')

##%######################################################%##
#                                                          #
####                   NMDS BACTERIA                    ####
#                                                          #
##%######################################################%##

library(vegan)
bacteria_dam_forest_otu_tab <- veganotu(bact_for_deh_opw)
bacteria_dam_forest_data <- data.frame(sample_data(bact_for_deh_opw))

# Set seed and make NMDS
set.seed(456782)
bacteria_dam_forest_NMDS <- metaMDS(bacteria_dam_forest_otu_tab, try = 500)
stressplot(bacteria_dam_forest_NMDS) #plot and check the stress of the NMDS ordination

# NMDS preparation for plotting
MDS1_bacteria_forest = bacteria_dam_forest_NMDS$points[,1]
MDS2_bacteria_forest = bacteria_dam_forest_NMDS$points[,2]
NMDS_bacteria_forest = data.frame(MDS1 = MDS1_bacteria_forest, MDS2 = MDS2_bacteria_forest, Manejo = bacteria_dam_forest_data$Manejo_Conteo)

# Reordering Land use factor variable
NMDS_bacteria_forest$Manejo <- factor(NMDS_bacteria_forest$Manejo, levels=c("Forest", "Open woodland", "Dehesa"))

# NDMS plot 
ggplot(NMDS_bacteria_forest, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(shape = Manejo), size = 5) +
  geom_vline(xintercept = 0, size = .1) +
  geom_hline(yintercept = 0, size = .1) +
  scale_shape_manual(values=c(21, 22, 24)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA)) +labs(title = "Bacterial core community") 



##%######################################################%##
#                                                          #
####                    NMDS HONGOS                     ####
#                                                          #
##%######################################################%##

library(vegan)
sort(sample_sums(mycota_for_deh_opw))
mycota_for_deh_opw <- prune_samples(sample_sums(mycota_for_deh_opw) >1000, mycota_for_deh_opw)
fungi_dam_forest_otu_tab <- veganotu(mycota_for_deh_opw)
fungi_dam_forest_data <- data.frame(sample_data(mycota_for_deh_opw))

# Set seed and make NMDS
set.seed(456782)
fungi_dam_forest_NMDS <- metaMDS(fungi_dam_forest_otu_tab, try = 500, k = 3)
stressplot(fungi_dam_forest_NMDS) #plot and check the stress of the NMDS ordination

# NMDS preparation for plotting
MDS1_fungi_forest = fungi_dam_forest_NMDS$points[,1]
MDS2_fungi_forest = fungi_dam_forest_NMDS$points[,2]
NMDS_fungi_forest = data.frame(MDS1 = MDS1_fungi_forest, MDS2 = MDS2_fungi_forest, Manejo = fungi_dam_forest_data$Manejo_Conteo)

## Reordering Land use factor variable
NMDS_fungi_forest$Manejo <- factor(NMDS_fungi_forest$Manejo, levels=c("Forest", "Open woodland", "Dehesa"))

# NDMS plot 
ggplot(NMDS_fungi_forest, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(shape = Manejo), size = 5) +
  geom_vline(xintercept = 0, size = .1) +
  geom_hline(yintercept = 0, size = .1) +
  scale_shape_manual(values=c(21, 22, 24)) +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(size = 0.2, 
                                 linetype = "solid"), plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = NA)) +labs(title = "Fungal core community") 
