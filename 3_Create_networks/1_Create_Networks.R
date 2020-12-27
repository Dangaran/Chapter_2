library(igraph)
library(SpiecEasi)
source('./accessory_functions/veganotu.R')
source('./accessory_functions/results_summary_network.R')


##%######################################################%##
#                                                          #
####                  Networks                          ####
#                                                          #
##%######################################################%##

# To calculate the networks we use the package SpiecEasi and to process the data igraph package
# https://github.com/zdk123/SpiecEasi

# Bacteria
# Forest
bact_network_for <- spiec.easi(illumina_forest_final, method = 'mb', lambda.min.ratio = 1e-2,
                               nlambda = 20, sel.criterion='bstars', pulsar.select = TRUE,
                               pulsar.params=list(rep.num = 50, ncores = 8, seed = 123))
# Dehesa
bact_network_deh <- spiec.easi(illumina_dehesa_final, method = 'mb', lambda.min.ratio = 1e-2,
                               nlambda = 20, sel.criterion='bstars', pulsar.select = TRUE,
                               pulsar.params=list(rep.num = 50, ncores = 8, seed = 123))
# Open woodland
bact_network_opw <- spiec.easi(illumina_openwood_final, method = 'mb', lambda.min.ratio = 1e-2,
                               nlambda = 20, sel.criterion='bstars', pulsar.select = TRUE,
                               pulsar.params=list(rep.num = 50, ncores = 8, seed = 123))
# Fungi
# Forest
fungi_network_for <- spiec.easi(mycota_forest_final, method = 'mb', lambda.min.ratio = 1e-2,
                                nlambda = 20, sel.criterion='bstars', pulsar.select = TRUE,
                                pulsar.params=list(rep.num = 50, ncores = 8, seed = 123))
# Dehesa
fungi_network_deh <- spiec.easi(mycota_dehesa_final, method = 'mb', lambda.min.ratio = 1e-2,
                                nlambda = 20, sel.criterion='bstars', pulsar.select = TRUE,
                                pulsar.params=list(rep.num = 50, ncores = 8, seed = 123))
# Open woodland
fungi_network_opw <- spiec.easi(mycota_openwood_final, method = 'mb', lambda.min.ratio = 1e-2,
                                nlambda = 20, sel.criterion='bstars', pulsar.select = TRUE,
                                pulsar.params=list(rep.num = 50, ncores = 8, seed = 123))



##%######################################################%##
#                                                          #
####               Create igraph objects                ####
#                                                          #
##%######################################################%##

# Create igraph object to analyze the networks
# Bacteria
bact_network_for <- unlist(bact_network_for)
bact_network_deh <- unlist(bact_network_deh)
bact_network_opw <- unlist(bact_network_opw)
# Fungi
fungi_network_for <- unlist(fungi_network_for)
fungi_network_deh <- unlist(fungi_network_deh)
fungi_network_opw <- unlist(fungi_network_opw)

# Add taxa_names to vertex
# Bacteria
bact_network_for_igraph <- adj2igraph(bact_network_for$refit, vertex.attr=list(name=taxa_names(illumina_forest_final)))
bact_network_deh_igraph <- adj2igraph(bact_network_deh$refit, vertex.attr=list(name=taxa_names(illumina_dehesa_final)))
bact_network_opw_igraph <- adj2igraph(bact_network_opw$refit, vertex.attr=list(name=taxa_names(illumina_openwood_final)))
# Fungi
fungi_network_for_igraph <- adj2igraph(fungi_network_for$refit, vertex.attr=list(name=taxa_names(mycota_forest_final)))
fungi_network_deh_igraph <- adj2igraph(fungi_network_deh$refit, vertex.attr=list(name=taxa_names(mycota_dehesa_final)))
fungi_network_opw_igraph <- adj2igraph(fungi_network_opw$refit, vertex.attr=list(name=taxa_names(mycota_openwood_final)))


# Create label object from the tax table in phyloseq objects
# Bacteria
bact_for_label <- data.frame(tax_table(illumina_forest_final)[,c(1:6)])
bact_deh_label <- data.frame(tax_table(illumina_dehesa_final)[,c(1:6)])
bact_opw_label <- data.frame(tax_table(illumina_openwood_final)[,c(1:6)])
# Fungi
fungi_for_label <- data.frame(tax_table(mycota_forest_final)[,c(1:6)])
fungi_deh_label <- data.frame(tax_table(mycota_dehesa_final)[,c(1:6)])
fungi_opw_label <- data.frame(tax_table(mycota_openwood_final)[,c(1:6)])


# Put label to nodes/vertex
# Bacteria
# Forest
V(bact_network_for_igraph)$Phylum <- as.vector(bact_for_label$Phylum)
V(bact_network_for_igraph)$Class <- as.vector(bact_for_label$Class)
V(bact_network_for_igraph)$Order <- as.vector(bact_for_label$Order)
V(bact_network_for_igraph)$Family <- as.vector(bact_for_label$Family)
V(bact_network_for_igraph)$Genus <- as.vector(bact_for_label$Genus)
V(bact_network_for_igraph)$OTU_ID <- as.vector(rownames(bact_for_label))
# Dehesa
V(bact_network_deh_igraph)$Phylum <- as.vector(bact_deh_label$Phylum)
V(bact_network_deh_igraph)$Class <- as.vector(bact_deh_label$Class)
V(bact_network_deh_igraph)$Order <- as.vector(bact_deh_label$Order)
V(bact_network_deh_igraph)$Family <- as.vector(bact_deh_label$Family)
V(bact_network_deh_igraph)$Genus <- as.vector(bact_deh_label$Genus)
V(bact_network_deh_igraph)$OTU_ID <- as.vector(rownames(bact_deh_label))
# Open woodland
V(bact_network_opw_igraph)$Phylum <- as.vector(bact_opw_label$Phylum)
V(bact_network_opw_igraph)$Class <- as.vector(bact_opw_label$Class)
V(bact_network_opw_igraph)$Order <- as.vector(bact_opw_label$Order)
V(bact_network_opw_igraph)$Family <- as.vector(bact_opw_label$Family)
V(bact_network_opw_igraph)$Genus <- as.vector(bact_opw_label$Genus)
V(bact_network_opw_igraph)$OTU_ID <- as.vector(rownames(bact_opw_label))

# Fungi
# Forest
V(fungi_network_for_igraph)$Phylum <- as.vector(fungi_for_label$Phylum)
V(fungi_network_for_igraph)$Class <- as.vector(fungi_for_label$Class)
V(fungi_network_for_igraph)$Order <- as.vector(fungi_for_label$Order)
V(fungi_network_for_igraph)$Family <- as.vector(fungi_for_label$Family)
V(fungi_network_for_igraph)$OTU_ID <- as.vector(rownames(fungi_for_label))
# Dehesa
V(fungi_network_deh_igraph)$Phylum <- as.vector(fungi_deh_label$Phylum)
V(fungi_network_deh_igraph)$Class <- as.vector(fungi_deh_label$Class)
V(fungi_network_deh_igraph)$Order <- as.vector(fungi_deh_label$Order)
V(fungi_network_deh_igraph)$Family <- as.vector(fungi_deh_label$Family)
V(fungi_network_deh_igraph)$OTU_ID <- as.vector(rownames(fungi_deh_label))
# Open woodland
V(fungi_network_opw_igraph)$Phylum <- as.vector(fungi_opw_label$Phylum)
V(fungi_network_opw_igraph)$Class <- as.vector(fungi_opw_label$Class)
V(fungi_network_opw_igraph)$Order <- as.vector(fungi_opw_label$Order)
V(fungi_network_opw_igraph)$Family <- as.vector(fungi_opw_label$Family)
V(fungi_network_opw_igraph)$OTU_ID <- as.vector(rownames(fungi_opw_label))


# Add Guild column according to OTU_ID name from the otu table for fungal communities
# Forest
Guild_for <- data.frame(Guild=funguild[match(names(V(fungi_network_for_igraph)), funguild$OTU.ID), 3])
V(fungi_network_for_igraph)$Guild <- as.vector(Guild_for$Guild)
# Dehesa
Guild_deh <- data.frame(Guild=funguild[match(names(V(fungi_network_deh_igraph)), funguild$OTU.ID), 3])
V(fungi_network_deh_igraph)$Guild <- as.vector(Guild_deh$Guild)
# Open woodland
Guild_opw <- data.frame(Guild=funguild[match(names(V(fungi_network_opw_igraph)), funguild$OTU.ID), 3])
V(fungi_network_opw_igraph)$Guild <- as.vector(Guild_opw$Guild)



##%######################################################%##
#                                                          #
####                  Summary results                   ####
#                                                          #
##%######################################################%##

# Bacteria
bact_result_for <- results_summary_network(bact_network_for_igraph)
bact_result_deh <- results_summary_network(bact_network_deh_igraph)
bact_result_opw <- results_summary_network(bact_network_opw_igraph)
# Fungi
fungi_result_for <- results_summary_network(fungi_network_for_igraph)
fungi_result_deh <- results_summary_network(fungi_network_deh_igraph)
fungi_result_opw <- results_summary_network(fungi_network_opw_igraph)

# Create dataframe with these data
Results_bact_fungi <- data.frame(matrix(ncol = 10, nrow = 6))
colnames(Results_bact_fungi) <- c("Network", "N_OTUs", "Larquest_clique", "3_nodes_cliques", "4_nodes_clique", "5_nodes_clique","Diameter", "Avg_path_length",
                                     "Connectance", "Clustering_coefficient")

Results_bact_fungi[1,1] <- "Bacteria_Forest"
Results_bact_fungi[1,2:10] <- bact_result_for
Results_bact_fungi[2,1] <- "Bacteria_Dehesa"
Results_bact_fungi[2,2:10] <- bact_result_deh
Results_bact_fungi[3,1] <- "Bacteria_Openwodland"
Results_bact_fungi[3,2:10] <- bact_result_opw

Results_bact_fungi[4,1] <- "Fungi_Forest"
Results_bact_fungi[4,2:10] <- fungi_result_for
Results_bact_fungi[5,1] <- "Fungi_Dehesa"
Results_bact_fungi[5,2:10] <- fungi_result_deh
Results_bact_fungi[6,1] <- "Fungi_Openwoodland"
Results_bact_fungi[6,2:10] <- fungi_result_opw

Results_bact_fungi
