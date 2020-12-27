source('./accessory_functions/plot_anova_diversity.R')
source('./accessory_functions/st_error_NA.R')


##%######################################################%##
#                                                          #
####             Create object to analyze               ####
#                                                          #
##%######################################################%##

# Extract otu_table, sample_data and tax_table for each land use and merge it together to analyze alpha and beta diversity

# Bacteria
# Forest
illumina_forest_final
for_otu_tab <- otu_table(illumina_forest_final)
for_sam_tab <- sample_data(illumina_forest_final)
for_tax_tab <- tax_table(illumina_forest_final)
bact_for_merge <- phyloseq(for_otu_tab, for_sam_tab, for_tax_tab)
# Dehesa
illumina_dehesa_final
deh_otu_tab <- otu_table(illumina_dehesa_final)
deh_sam_tab <- sample_data(illumina_dehesa_final)
deh_tax_tab <- tax_table(illumina_dehesa_final)
bact_deh_merge <- phyloseq(deh_otu_tab, deh_sam_tab, deh_tax_tab)
# Open woodland
illumina_openwood_final
opw_otu_tab <- otu_table(illumina_openwood_final)
opw_sam_tab <- sample_data(illumina_openwood_final)
opw_tax_tab <- tax_table(illumina_openwood_final)
bact_opw_merge <- phyloseq(opw_otu_tab, opw_sam_tab, opw_tax_tab)

# Merge files
bact_for_deh <- merge_phyloseq(bact_for_merge, bact_deh_merge)
bact_for_deh_opw <- merge_phyloseq(bact_for_merge, bact_deh_merge,bact_opw_merge)
bact_for_deh_opw

# Fungi
# As fungal data doesn't have a taxonomic tree is not necessary to extract the data 
mycota_forest_final
mycota_dehesa_final
mycota_openwood_final
# Merge files
mycota_for_deh <- merge_phyloseq(mycota_forest_final, mycota_dehesa_final)
mycota_for_deh_opw <- merge_phyloseq(mycota_forest_final, mycota_dehesa_final,mycota_openwood_final)
mycota_for_deh_opw


##%######################################################%##
#                                                          #
####             Alpha diversity - Shannon              ####
#                                                          #
##%######################################################%##

# Bacteria
# Calculate index and extract data
shannon_bact <- plot_anova_diversity(bact_for_deh_opw, method = c("shannon"), grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
shannon_bact_data <- shannon_bact$data
# Anova analysis
anova_bact <- aov(value ~ Manejo_Conteo, shannon_bact_data)
summary(anova_bact)
TukeyHSD(anova_bact)
# Calculate mean and standard error
aggregate(value ~ Manejo_Conteo, shannon_bact_data, mean)
aggregate(value ~ Manejo_Conteo, shannon_bact_data, st_error_NA)

# Fungi
# Calculate index and extract data
shannon_fungi <- plot_anova_diversity(mycota_for_deh_opw, method = c("shannon"), grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
shannon_fungi_data <- shannon_fungi$data
# Anova analysis
anova_fungi <- aov(value ~ Manejo_Conteo, shannon_fungi_data)
summary(anova_fungi)
TukeyHSD(anova_fungi)
# Calculate mean and standard error
aggregate(value ~ Manejo_Conteo, shannon_fungi_data, mean)
aggregate(value ~ Manejo_Conteo, shannon_fungi_data, st_error_NA)



##%######################################################%##
#                                                          #
####             Alpha diversity - Simpson              ####
#                                                          #
##%######################################################%##

# Bacteria
# Calculate index and extract data
simpson_bact <- plot_anova_diversity(bact_for_deh_opw, method = c("simpson"), grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
simpson_bact_data <- simpson_bact$data
# Anova analysis
anova_bact <- aov(value ~ Manejo_Conteo, simpson_bact_data)
summary(anova_bact)
TukeyHSD(anova_bact)
# Calculate mean and standard error
aggregate(value ~ Manejo_Conteo, simpson_bact_data, mean)
aggregate(value ~ Manejo_Conteo, simpson_bact_data, st_error_NA)

# Fungi
# Calculate index and extract data
simpson_fungi <- plot_anova_diversity(mycota_for_deh_opw, method = c("simpson"), grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
simpson_fungi_data <- simpson_fungi$data
# Anova analysis
anova_fungi <- aov(value ~ Manejo_Conteo, simpson_fungi_data)
summary(anova_fungi)
TukeyHSD(anova_fungi)
# Calculate mean and standard error
aggregate(value ~ Manejo_Conteo, simpson_fungi_data, mean)
aggregate(value ~ Manejo_Conteo, simpson_fungi_data, st_error_NA)



##%######################################################%##
#                                                          #
####        Alpha diversity - Pielou’s evenness         ####
#                                                          #
##%######################################################%##

# Bacteria
# Calculate index and extract data
evenness_bact <- plot_anova_diversity(bact_for_deh_opw, method = c("evenness"), grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
evenness_bact_data <- evenness_bact$data
# Anova analysis
anova_bact <- aov(value ~ Manejo_Conteo, evenness_bact_data)
summary(anova_bact)
TukeyHSD(anova_bact)
# Calculate mean and standard error
aggregate(value ~ Manejo_Conteo, evenness_bact_data, mean)
aggregate(value ~ Manejo_Conteo, evenness_bact_data, st_error_NA)

# Fungi
# Calculate index and extract data
evenness_fungi <- plot_anova_diversity(mycota_for_deh_opw, method = c("evenness"), grouping_column = "Manejo_Conteo", pValueCutoff = 0.05)
evenness_fungi_data <- evenness_fungi$data
# Anova analysis
anova_fungi <- aov(value ~ Manejo_Conteo, evenness_fungi_data)
summary(anova_fungi)
TukeyHSD(anova_fungi)
# Calculate mean and standard error
aggregate(value ~ Manejo_Conteo, evenness_fungi_data, mean)
aggregate(value ~ Manejo_Conteo, evenness_fungi_data, st_error_NA)