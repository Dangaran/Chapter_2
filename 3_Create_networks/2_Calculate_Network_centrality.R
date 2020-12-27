library(igraph)

##%######################################################%##
#                                                          #
####               Betweenness Centrality               ####
#                                                          #
##%######################################################%##

# Find the hub OTUs
# Compute node degrees (#links) and use that to set node size

# Bacteria
# Forest
centr_betw_illumina_for <- centr_betw(bact_network_for_igraph)$res
V(bact_network_for_igraph)$betw <- centr_betw_illumina_for/max(centr_betw_illumina_for)
centr_betw_illumina_for <- data.frame(ID = V(bact_network_for_igraph)$OTU_ID, centr_betw = V(bact_network_for_igraph)$betw, Phylum = V(bact_network_for_igraph)$Phylum, 
                                      Class = V(bact_network_for_igraph)$Class, Order = V(bact_network_for_igraph)$Order,
                                      Family = V(bact_network_for_igraph)$Family)
centr_betw_illumina_for <- centr_betw_illumina_for[order(-centr_betw_illumina_for$centr_betw),]
head(centr_betw_illumina_for,15)
# Dehesa
centr_betw_illumina_deh <- centr_betw(bact_network_deh_igraph)$res
V(bact_network_deh_igraph)$betw <- centr_betw_illumina_deh/max(centr_betw_illumina_deh)
centr_betw_illumina_deh <- data.frame(ID = V(bact_network_deh_igraph)$OTU_ID, centr_betw = V(bact_network_deh_igraph)$betw, Phylum = V(bact_network_deh_igraph)$Phylum, 
                                      Class = V(bact_network_deh_igraph)$Class, Order = V(bact_network_deh_igraph)$Order,
                                      Family = V(bact_network_deh_igraph)$Family)
centr_betw_illumina_deh <- centr_betw_illumina_deh[order(-centr_betw_illumina_deh$centr_betw),]
head(centr_betw_illumina_deh,15)
# Open woodland
centr_betw_illumina_opw <- centr_betw(bact_network_opw_igraph)$res
V(bact_network_opw_igraph)$betw <- centr_betw_illumina_opw/max(centr_betw_illumina_opw)
centr_betw_illumina_opw <- data.frame(ID = V(bact_network_opw_igraph)$OTU_ID, centr_betw = V(bact_network_opw_igraph)$betw, Phylum = V(bact_network_opw_igraph)$Phylum, 
                                      Class = V(bact_network_opw_igraph)$Class, Order = V(bact_network_opw_igraph)$Order,
                                      Family = V(bact_network_opw_igraph)$Family)
centr_betw_illumina_opw <- centr_betw_illumina_opw[order(-centr_betw_illumina_opw$centr_betw),]
head(centr_betw_illumina_opw,15)

# Fungi
# Forest
centr_betw_mycota_for <- centr_betw(fungi_network_for_igraph)$res
V(fungi_network_for_igraph)$betw <- centr_betw_mycota_for/max(centr_betw_mycota_for)
centr_betw_mycota_for <- data.frame(ID = V(fungi_network_for_igraph)$OTU_ID, centr_betw = V(fungi_network_for_igraph)$betw, Phylum = V(fungi_network_for_igraph)$Phylum, 
                                    Class = V(fungi_network_for_igraph)$Class, Order = V(fungi_network_for_igraph)$Order,
                                    Family = V(fungi_network_for_igraph)$Family, Guild = V(fungi_network_for_igraph)$Guild)
centr_betw_mycota_for <- centr_betw_mycota_for[order(-centr_betw_mycota_for$centr_betw),]
head(centr_betw_mycota_for,15)
# Dehesa
centr_betw_mycota_deh <- centr_betw(fungi_network_deh_igraph)$res
V(fungi_network_deh_igraph)$betw <- centr_betw_mycota_deh/max(centr_betw_mycota_deh)
centr_betw_mycota_deh <- data.frame(ID = V(fungi_network_deh_igraph)$OTU_ID, centr_betw = V(fungi_network_deh_igraph)$betw, Phylum = V(fungi_network_deh_igraph)$Phylum, 
                                    Class = V(fungi_network_deh_igraph)$Class, Order = V(fungi_network_deh_igraph)$Order,
                                    Family = V(fungi_network_deh_igraph)$Family, Guild = V(fungi_network_deh_igraph)$Guild)
centr_betw_mycota_deh <- centr_betw_mycota_deh[order(-centr_betw_mycota_deh$centr_betw),]
head(centr_betw_mycota_deh,15)
# Open woodland
centr_betw_mycota_opw <- centr_betw(fungi_network_opw_igraph)$res
V(fungi_network_opw_igraph)$betw <- centr_betw_mycota_opw/max(centr_betw_mycota_opw)
centr_betw_mycota_opw <- data.frame(ID = V(fungi_network_opw_igraph)$OTU_ID, centr_betw = V(fungi_network_opw_igraph)$betw, Phylum = V(fungi_network_opw_igraph)$Phylum, 
                                    Class = V(fungi_network_opw_igraph)$Class, Order = V(fungi_network_opw_igraph)$Order,
                                    Family = V(fungi_network_opw_igraph)$Family, Guild = V(fungi_network_opw_igraph)$Guild)
centr_betw_mycota_opw <- centr_betw_mycota_opw[order(-centr_betw_mycota_opw$centr_betw),]
head(centr_betw_mycota_opw,15)


##%######################################################%##
#                                                          #
####                 Degree Centrality                  ####
#                                                          #
##%######################################################%##

# Bacteria
# Forest
degree_illumina_for <- data.frame(ID = V(bact_network_for_igraph)$OTU_ID, degree = degree(bact_network_for_igraph)/max(degree(bact_network_for_igraph)),
                                  Phylum = V(bact_network_for_igraph)$Phylum, 
                                  Class = V(bact_network_for_igraph)$Class, 
                                  Order = V(bact_network_for_igraph)$Order,
                                  Family = V(bact_network_for_igraph)$Family) 
degree_illumina_for <- degree_illumina_for[order(-degree_illumina_for$degree),]
head(degree_illumina_for)
# Dehesa
degree_illumina_deh <- data.frame(ID = V(bact_network_deh_igraph)$OTU_ID, degree = degree(bact_network_deh_igraph)/max(degree(bact_network_deh_igraph)),
                                  Phylum = V(bact_network_deh_igraph)$Phylum, 
                                  Class = V(bact_network_deh_igraph)$Class, 
                                  Order = V(bact_network_deh_igraph)$Order,
                                  Family = V(bact_network_deh_igraph)$Family) 
degree_illumina_deh <- degree_illumina_deh[order(-degree_illumina_deh$degree),]
degree_illumina_deh
# Open woodland
degree_illumina_opw <- data.frame(ID = V(bact_network_opw_igraph)$OTU_ID, degree = degree(bact_network_opw_igraph)/max(degree(bact_network_opw_igraph)),
                                  Phylum = V(bact_network_opw_igraph)$Phylum, 
                                  Class = V(bact_network_opw_igraph)$Class, 
                                  Order = V(bact_network_opw_igraph)$Order,
                                  Family = V(bact_network_opw_igraph)$Family) 
degree_illumina_opw <- degree_illumina_opw[order(-degree_illumina_opw$degree),]
degree_illumina_opw

# Fungal
# Forest
degree_mycota_for <- data.frame(ID = V(fungi_network_for_igraph)$OTU_ID, degree = degree(fungi_network_for_igraph)/max(degree(fungi_network_for_igraph)),
                                Phylum = V(fungi_network_for_igraph)$Phylum, 
                                Class = V(fungi_network_for_igraph)$Class, 
                                Order = V(fungi_network_for_igraph)$Order,
                                Family = V(fungi_network_for_igraph)$Family,
                                Guild = V(fungi_network_for_igraph)$Guild) 
degree_mycota_for <- degree_mycota_for[order(-degree_mycota_for$degree),]
degree_mycota_for
# Dehesa
degree_mycota_deh <- data.frame(ID = V(fungi_network_deh_igraph)$OTU_ID, degree = degree(fungi_network_deh_igraph)/max(degree(fungi_network_deh_igraph)),
                                Phylum = V(fungi_network_deh_igraph)$Phylum, 
                                Class = V(fungi_network_deh_igraph)$Class, 
                                Order = V(fungi_network_deh_igraph)$Order,
                                Family = V(fungi_network_deh_igraph)$Family,
                                Guild = V(fungi_network_deh_igraph)$Guild) 
degree_mycota_deh <- degree_mycota_deh[order(-degree_mycota_deh$degree),]
degree_mycota_deh
# Open woodland
degree_mycota_opw <- data.frame(ID = V(fungi_network_opw_igraph)$OTU_ID, degree = degree(fungi_network_opw_igraph)/max(degree(fungi_network_opw_igraph)),
                                Phylum = V(fungi_network_opw_igraph)$Phylum, 
                                Class = V(fungi_network_opw_igraph)$Class, 
                                Order = V(fungi_network_opw_igraph)$Order,
                                Family = V(fungi_network_opw_igraph)$Family,
                                Guild = V(fungi_network_opw_igraph)$Guild) 
degree_mycota_opw <- degree_mycota_opw[order(-degree_mycota_opw$degree),]
degree_mycota_opw

##%######################################################%##
#                                                          #
####               Closeness Centrality                 ####
#                                                          #
##%######################################################%##

# Bacteria
# Forest
centr_clo_illumina_for <- centr_clo(bact_network_for_igraph)$res
centr_clo_illumina_for <- data.frame(ID = V(bact_network_for_igraph)$OTU_ID, centr_clo = centr_clo_illumina_for, Phylum = V(bact_network_for_igraph)$Phylum, 
                                     Class = V(bact_network_for_igraph)$Class, Order = V(bact_network_for_igraph)$Order,
                                     Family = V(bact_network_for_igraph)$Family)
centr_clo_illumina_for <- centr_clo_illumina_for[order(-centr_clo_illumina_for$centr_clo),]
head(centr_clo_illumina_for,15)
# Dehesa
centr_clo_illumina_deh <- centr_clo(bact_network_deh_igraph)$res
centr_clo_illumina_deh <- data.frame(ID = V(bact_network_deh_igraph)$OTU_ID, centr_clo = centr_clo_illumina_deh, Phylum = V(bact_network_deh_igraph)$Phylum, 
                                     Class = V(bact_network_deh_igraph)$Class, Order = V(bact_network_deh_igraph)$Order,
                                     Family = V(bact_network_deh_igraph)$Family)
centr_clo_illumina_deh <- centr_clo_illumina_deh[order(-centr_clo_illumina_deh$centr_clo),]
head(centr_clo_illumina_deh,15)
# Open woodland
centr_clo_illumina_opw <- centr_clo(bact_network_opw_igraph)$res
centr_clo_illumina_opw <- data.frame(ID = V(bact_network_opw_igraph)$OTU_ID, centr_clo = centr_clo_illumina_opw, Phylum = V(bact_network_opw_igraph)$Phylum, 
                                     Class = V(bact_network_opw_igraph)$Class, Order = V(bact_network_opw_igraph)$Order,
                                     Family = V(bact_network_opw_igraph)$Family)
centr_clo_illumina_opw <- centr_clo_illumina_opw[order(-centr_clo_illumina_opw$centr_clo),]
head(centr_clo_illumina_opw,15)

# Fungal
# Forest
centr_clo_mycota_for <- centr_clo(fungi_network_for_igraph)$res
centr_clo_mycota_for <- data.frame(ID = V(fungi_network_for_igraph)$OTU_ID, centr_clo = centr_clo_mycota_for, Phylum = V(fungi_network_for_igraph)$Phylum, 
                                   Class = V(fungi_network_for_igraph)$Class, Order = V(fungi_network_for_igraph)$Order,
                                   Family = V(fungi_network_for_igraph)$Family, Guild = V(fungi_network_for_igraph)$Guild)
centr_clo_mycota_for <- centr_clo_mycota_for[order(-centr_clo_mycota_for$centr_clo),]
head(centr_clo_mycota_for,15)
# Dehesa
centr_clo_mycota_deh <- centr_clo(fungi_network_deh_igraph)$res
centr_clo_mycota_deh <- data.frame(ID = V(fungi_network_deh_igraph)$OTU_ID, centr_clo = centr_clo_mycota_deh, Phylum = V(fungi_network_deh_igraph)$Phylum, 
                                   Class = V(fungi_network_deh_igraph)$Class, Order = V(fungi_network_deh_igraph)$Order,
                                   Family = V(fungi_network_deh_igraph)$Family, Guild = V(fungi_network_deh_igraph)$Guild)
centr_clo_mycota_deh <- centr_clo_mycota_deh[order(-centr_clo_mycota_deh$centr_clo),]
head(centr_clo_mycota_deh,15)
# Open woodland
centr_clo_mycota_opw <- centr_clo(fungi_network_opw_igraph)$res
centr_clo_mycota_opw <- data.frame(ID = V(fungi_network_opw_igraph)$OTU_ID, centr_clo = centr_clo_mycota_opw, Phylum = V(fungi_network_opw_igraph)$Phylum, 
                                   Class = V(fungi_network_opw_igraph)$Class, Order = V(fungi_network_opw_igraph)$Order,
                                   Family = V(fungi_network_opw_igraph)$Family, Guild = V(fungi_network_opw_igraph)$Guild)
centr_clo_mycota_opw <- centr_clo_mycota_opw[order(-centr_clo_mycota_opw$centr_clo),]
head(centr_clo_mycota_opw,15)


##%######################################################%##
#                                                          #
####                    Transitivity                    ####
#                                                          #
##%######################################################%##

# Bacteria
# Forest
trans_illumina_for <- data.frame(ID = V(bact_network_for_igraph)$OTU_ID, trans = transitivity(bact_network_for_igraph, type = "local"))
trans_illumina_for <- trans_illumina_for[order(-trans_illumina_for$trans),]
# Dehesa
trans_illumina_deh <- data.frame(ID = V(bact_network_deh_igraph)$OTU_ID, trans = transitivity(bact_network_deh_igraph, type = "local"))
trans_illumina_deh <- trans_illumina_deh[order(-trans_illumina_deh$trans),]
# Open woodland
trans_illumina_opw <- data.frame(ID = V(bact_network_opw_igraph)$OTU_ID, trans = transitivity(bact_network_opw_igraph, type = "local"))
trans_illumina_opw <- trans_illumina_opw[order(-trans_illumina_opw$trans),]

# Fungi
# Forest
trans_mycota_for <- data.frame(ID = V(fungi_network_for_igraph)$OTU_ID, trans = transitivity(fungi_network_for_igraph, type = "local"))
trans_mycota_for <- trans_mycota_for[order(-trans_mycota_for$trans),]
# Dehesa
trans_mycota_deh <- data.frame(ID = V(fungi_network_deh_igraph)$OTU_ID, trans = transitivity(fungi_network_deh_igraph, type = "local"))
trans_mycota_deh <- trans_mycota_deh[order(-trans_mycota_deh$trans),]
# Open woodland
trans_mycota_opw <- data.frame(ID = V(fungi_network_opw_igraph)$OTU_ID, trans = transitivity(fungi_network_opw_igraph, type = "local"))
trans_mycota_opw <- trans_mycota_opw[order(-trans_mycota_opw$trans),]
