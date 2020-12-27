library(phyloseq)
library(ggplot2)
source('./accessory_functions/remove_bad_coverage_samples.R')

##%######################################################%##
#                                                          #
####                Recommended websites                ####
#                                                          #
##%######################################################%##

# http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
# https://github.com/microbiome/microbiome/tree/master/R
# http://genoweb.toulouse.inra.fr/~formation/15_FROGS/5-June2016/FROGS_phyloseq_23062016.pdf
# http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html


##%######################################################%##
#                                                          #
####                Prepare input file                  ####
#                                                          #
##%######################################################%##

illumina_fungi <- import_biom("Datos/Fungi/otu_table.json", parseFunction = parse_taxonomy_default) 
colnames(tax_table(illumina_fungi)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Specie")

# Add sample variables
Sample.data <- read.csv(file = "Datos/Bacteria/base_datos_metro.csv", row.names = 1, header = TRUE, sep = ";")
Sample.data$Manejo_Conteo <- factor(Sample.data$Manejo_Conteo, levels = c("Dehesa", "Open woodland", "Forest"))
Sample.data$Type <- factor(Sample.data$Type, levels = c("Healthy", "Affected", "Dead", "Bare soil"))
sample_data(illumina_fungi) <- Sample.data


##%######################################################%##
#                                                          #
####                 Sequencing filters                 ####
#                                                          #
##%######################################################%##

# Good's coverage: Singletons/total reads (First filter)
remove_seq <- remove.bad.coverage.samples(illumina_fungi)
remove_seq$coverage_plot

illumina_fungi_good <- remove.bad.coverage.samples(illumina_fungi,85)

# Remove taxa with less than 10 reads (Second filter) http://fungal-sequencing-methods-discussion.blogspot.dk/p/funguild.html
illumina_fungi_10seq <-  prune_taxa(taxa_sums(illumina_fungi_good) > 10, illumina_fungi_good)


# Remove non-fungal taxa
mycota <- subset_taxa(illumina_fungi_10seq, Phylum=="p__Basidiomycota" | Phylum=="p__Ascomycota" | Phylum=="p__Chytridiomycota" | Phylum=="p__Mucoromycota"
                      | Phylum=="p__Mortierellomycota" | Phylum=="p__GS19" | Phylum=="p__Rozellomycota" | Phylum=="p__Entomophthoromycota"
                      | Phylum=="p__Monoblepharomycota" | Phylum=="p__Calcarisporiellomycota" | Phylum=="p__Zoopagomycota" | Phylum=="p__Glomeromycota"
                      | Phylum=="p__Oomycota" | Phylum=="p__Kickxellomycota" | Phylum=="p__Olpidiomycota" | Phylum=="p__Neocallimastigomycota"
                      | Phylum=="p__Entorrhizomycota" | Phylum=="p__Blastocladiomycota" | Phylum=="p__Aphelidiomycota" | Phylum=="p__Bacillariophyta")

hist(sample_sums(mycota))
sum(sample_sums(mycota))
sort(sample_sums(mycota))

# Remove samples with less than 1000 sequences
mycota <- prune_samples(sample_sums(mycota)>2000, mycota)

# Rename species s__unidentified 
tax_table(mycota)@.Data[,7] <- ifelse(tax_table(mycota)@.Data[,7] == "s__unidentified", sprintf("s__unidentified%d", 1:6625), tax_table(mycota)@.Data[,7])

mycota_glom <- tax_glom(mycota, taxrank = "Specie")


##%######################################################%##
#                                                          #
####      Info about physeq data after filtering        ####
#                                                          #
##%######################################################%##

# Number of sequences per sample
sort(sample_sums(mycota_glom))
hist(sample_sums(mycota_glom))

sort(sample_sums(mycota_glom))

# Mean, median, minimum and maximum reads per sample
summary(sample_sums(mycota_glom))


##%######################################################%##
#                                                          #
####       Remove OTUs with a relative abundance        ####
####          higher than 0.01 in all samples           ####
#                                                          #
##%######################################################%##

fungi_damaged <- subset_samples(mycota_glom, Treatment == "Damaged") # Remove control samples

minTotRelAbun = 0.0001
x_fun = taxa_sums(fungi_damaged)
keepTaxa_fun = which((x_fun / sum(x_fun)) > minTotRelAbun)
fun_damaged_relab = prune_taxa(names(keepTaxa_fun), fungi_damaged)
fun_damaged_relab
