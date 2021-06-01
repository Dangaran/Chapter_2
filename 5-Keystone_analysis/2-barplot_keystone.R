# Bacteria
database_keyst_bact <- read.csv("Datos/database_keystone_phyla.csv", h = T, sep = ";")
## Reordering database_keyst_bact$Land_use
database_keyst_bact$Land_use <- factor(database_keyst_bact$Land_use, levels=c("Forest", "Open woodland", "Dehesa"))

library(reshape2)
library(RColorBrewer)
database_keyst_bact <- melt(database_keyst_bact)
database_keyst_bact_for <- subset(database_keyst_bact, Land_use == "Forest")
database_keyst_bact_opw <- subset(database_keyst_bact, Land_use == "Open woodland")
database_keyst_bact_deh <- subset(database_keyst_bact, Land_use == "Dehesa")

database_keyst_bact_for$Percentage <- (database_keyst_bact_for$value/sum(database_keyst_bact_for$value))*100
database_keyst_bact_opw$Percentage <- (database_keyst_bact_opw$value/sum(database_keyst_bact_opw$value))*100
database_keyst_bact_deh$Percentage <- (database_keyst_bact_deh$value/sum(database_keyst_bact_deh$value))*100

database_keyst_bact <- rbind(database_keyst_bact_for, database_keyst_bact_opw, database_keyst_bact_deh)
library(ggplot2)
ggplot(database_keyst_bact, aes(variable, Percentage,  fill = Phylum)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(brewer.pal(8,"Set2"))) +
  facet_wrap(~Land_use) +
  labs(x = "Keystone type", y = "Percentage (%)") + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1)) + theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                                                            plot.caption = element_text(vjust = 1),
                                                                                                            text=element_text(family="Times New Roman", face="plain", size=19),
                                                                                                            panel.background = element_rect(fill = "gray95"))


# Fungi
database_keyst_fungi <- read.csv("Datos/database_keystone_phyla_fungi.csv", h = T, sep = ";")
## Reordering database_keyst_fungi$Land_use
database_keyst_fungi$Land_use <- factor(database_keyst_fungi$Land_use, levels=c("Forest", "Open woodland", "Dehesa"))

library(reshape2)
database_keyst_fungi <- melt(database_keyst_fungi)
database_keyst_fungi_for <- subset(database_keyst_fungi, Land_use == "Forest")
database_keyst_fungi_opw <- subset(database_keyst_fungi, Land_use == "Open woodland")
database_keyst_fungi_deh <- subset(database_keyst_fungi, Land_use == "Dehesa")

database_keyst_fungi_for$Percentage <- (database_keyst_fungi_for$value/sum(database_keyst_fungi_for$value))*100
database_keyst_fungi_opw$Percentage <- (database_keyst_fungi_opw$value/sum(database_keyst_fungi_opw$value))*100
database_keyst_fungi_deh$Percentage <- (database_keyst_fungi_deh$value/sum(database_keyst_fungi_deh$value))*100

database_keyst_fungi <- rbind(database_keyst_fungi_for, database_keyst_fungi_opw, database_keyst_fungi_deh)

library(ggplot2)
ggplot(database_keyst_fungi, aes(variable, Percentage,  fill = Phylum)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(brewer.pal(4,"Set2"))) +
  facet_wrap(~Land_use) +
  labs(x = "Keystone type", y = "Percentage (%)") + theme(plot.subtitle = element_text(vjust = 1), 
                                                   plot.caption = element_text(vjust = 1))+ 
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),text=element_text(family="Times New Roman", face="plain", size=19),
                                                                 panel.background = element_rect(fill = "gray95"))

