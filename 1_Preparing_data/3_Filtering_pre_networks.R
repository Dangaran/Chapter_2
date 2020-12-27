##### Data from the last step in previous scripts


##%######################################################%##
#                                                          #
####        Calculate number of reads per OTUs          ####
#                                                          #
##%######################################################%##

# Forest
# Bacteria
illumina_forest <- subset_samples(bact_damaged_relab, Manejo_Conteo == "Forest")
mean(taxa_sums(illumina_forest))
# Fungi
mycota_forest <- subset_samples(fun_damaged_relab, Manejo_Conteo == "Forest")
mean(taxa_sums(mycota_forest))

# Dehesa
# Bacteria
illumina_dehesa <- subset_samples(bact_damaged_relab, Manejo_Conteo == "Dehesa")
mean(taxa_sums(illumina_dehesa))
# Fungi
mycota_dehesa <- subset_samples(fun_damaged_relab, Manejo_Conteo == "Dehesa")
mean(taxa_sums(mycota_dehesa))

# Open woodland
# Bacteria
illumina_openwood <- subset_samples(bact_damaged_relab, Manejo_Conteo == "Open woodland")
mean(taxa_sums(illumina_openwood))
# Fungi
mycota_openwood <- subset_samples(fun_damaged_relab, Manejo_Conteo == "Open woodland")
mean(taxa_sums(mycota_openwood))


##%######################################################%##
#                                                          #
####           Filter taxa per number of read           ####
#                                                          #
##%######################################################%##

# First filter: OTUs must have at least 400 reads for bacteria and 200 reads for fungi to be considered as core community

# Bacteria
illumina_forest.filter <- prune_taxa(taxa_sums(illumina_forest) > 400, illumina_forest)
illumina_dehesa.filter <- prune_taxa(taxa_sums(illumina_dehesa) > 400, illumina_dehesa)
illumina_openwood.filter <- prune_taxa(taxa_sums(illumina_openwood) > 400, illumina_openwood)

# Fungi
mycota_forest.filter <- prune_taxa(taxa_sums(mycota_forest) > 200, mycota_forest)
mycota_dehesa.filter <- prune_taxa(taxa_sums(mycota_dehesa) > 200, mycota_dehesa)
mycota_openwood.filter <- prune_taxa(taxa_sums(mycota_openwood) > 200, mycota_openwood)


##%######################################################%##
#                                                          #
####                Filter taxa per site                ####
#                                                          #
##%######################################################%##

# Second filter: OTUs must appear in at least 3 different Sites to be considered as core community

# Bacteria
# Forest
illumina_forest.site <- merge_samples(illumina_forest.filter, group = "Site")
illumina_forest.otus <- filter_taxa(illumina_forest.site, function(x) sum(x >= 1) >= (3), TRUE)
illumina_forest.otus
otu_names_illum_for <- taxa_names(illumina_forest.otus)
illumina_forest_final <- prune_taxa(otu_names_illum_for, illumina_forest.filter)
illumina_forest_final
# Dehesa
illumina_dehesa.site <- merge_samples(illumina_dehesa.filter, group = "Site")
illumina_dehesa.otus <- filter_taxa(illumina_dehesa.site, function(x) sum(x >= 1) >= (3), TRUE)
illumina_dehesa.otus
otu_names_illum_for <- taxa_names(illumina_dehesa.otus)
illumina_dehesa_final <- prune_taxa(otu_names_illum_for, illumina_dehesa.filter)
illumina_dehesa_final
# Open woodland
illumina_openwood.site <- merge_samples(illumina_openwood.filter, group = "Site")
illumina_openwood.otus <- filter_taxa(illumina_openwood.site, function(x) sum(x >= 1) >= (3), TRUE)
illumina_openwood.otus
otu_names_illum_for <- taxa_names(illumina_openwood.otus)
illumina_openwood_final <- prune_taxa(otu_names_illum_for, illumina_openwood.filter)
illumina_openwood_final


# Fungi
# Forest
mycota_forest.site <- merge_samples(mycota_forest.filter, group = "Site")
mycota_forest.otus <- filter_taxa(mycota_forest.site, function(x) sum(x >= 1) >= (3), TRUE)
mycota_forest.otus
otu_names_illum_for <- taxa_names(mycota_forest.otus)
mycota_forest_final <- prune_taxa(otu_names_illum_for, mycota_forest.filter)
mycota_forest_final
# Dehesa
mycota_dehesa.site <- merge_samples(mycota_dehesa.filter, group = "Site")
mycota_dehesa.otus <- filter_taxa(mycota_dehesa.site, function(x) sum(x >= 1) >= (3), TRUE)
mycota_dehesa.otus
otu_names_illum_for <- taxa_names(mycota_dehesa.otus)
mycota_dehesa_final <- prune_taxa(otu_names_illum_for, mycota_dehesa.filter)
mycota_dehesa_final
# Open woodland
mycota_openwood.site <- merge_samples(mycota_openwood.filter, group = "Site")
mycota_openwood.otus <- filter_taxa(mycota_openwood.site, function(x) sum(x >= 1) >= (3), TRUE)
mycota_openwood.otus
otu_names_illum_for <- taxa_names(mycota_openwood.otus)
mycota_openwood_final <- prune_taxa(otu_names_illum_for, mycota_openwood.filter)
mycota_openwood_final

