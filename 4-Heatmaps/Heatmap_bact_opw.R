relab_bact_opw_core <- transform_sample_counts(illumina_openwood_final, function(x) x / sum(x))
relab_bact_opw_core

# Merge by phylum filo
relab_bact_opw_core <-  tax_glom(relab_bact_opw_core, "Phylum")

# Add tree influence as variable
core_bact_opw_data <- data.frame(sample_data(relab_bact_opw_core))
core_bact_opw_data$Tree_inf <- core_bact_opw_data$Porcentaje_defol
core_bact_opw_data$Tree_inf <- ifelse(is.na(core_bact_opw_data$Tree_inf) == TRUE, 150, core_bact_opw_data$Tree_inf)
colnames(core_bact_opw_data)[c(15,30,31, 36, 43,46,45,50,49,53)] <- c("Defol", "Temp", "Prec","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH")
sample_data(relab_bact_opw_core) <- core_bact_opw_data


# MicrobiomeSEQ -----
taxa_names(relab_bact_opw_core) <- paste("OTU_ID", seq(1,length(as.vector(tax_table(relab_bact_opw_core)[,2]))), sep = "_")
colnames(core_bact_opw_data)
correlacion_nut_func <- taxa.env.correlation(relab_bact_opw_core, grouping_column = "Manejo_Conteo", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                             select.variables = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"))

# Add taxonomy to df
relab_bact_opw_core_tax_table <- tax_table(relab_bact_opw_core)
correlacion_nut_func$Phylum <- as.factor(relab_bact_opw_core_tax_table[match(correlacion_nut_func$Taxa, rownames(relab_bact_opw_core_tax_table)), 2])

# Sort variables
correlacion_nut_func$Env <- factor(correlacion_nut_func$Env, levels=c("Temp", "Prec", "pH", "Defol", "Tree_inf", "Org_C", "NH4_N", "NO3_N", "Av_P", "RH", "Ramm", "Rnit"))

# Create dummy phyla
AD3 <- data.frame(Taxa = rep("Dummy",12), 
                  Env = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"),
                  Correlation = rep(0,12),
                  Pvalue = rep(0,12),
                  Type = rep("Open woodland",12),
                  AdjPvalue = rep(0,12),
                  Significance = rep(NA,12),
                  Phylum = rep("AD3",12))

Armatimonadetes <- data.frame(Taxa = rep("Dummy",12), 
                              Env = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"),
                              Correlation = rep(0,12),
                              Pvalue = rep(0,12),
                              Type = rep("Open woodland",12),
                              AdjPvalue = rep(0,12),
                              Significance = rep(NA,12),
                              Phylum = rep("Armatimonadetes",12))

Cyanobacteria <- data.frame(Taxa = rep("Dummy",12), 
                            Env = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"),
                            Correlation = rep(0,12),
                            Pvalue = rep(0,12),
                            Type = rep("Open woodland",12),
                            AdjPvalue = rep(0,12),
                            Significance = rep(NA,12),
                            Phylum = rep("Cyanobacteria",12))

Elusimicrobia <- data.frame(Taxa = rep("Dummy",12), 
                            Env = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"),
                            Correlation = rep(0,12),
                            Pvalue = rep(0,12),
                            Type = rep("Open woodland",12),
                            AdjPvalue = rep(0,12),
                            Significance = rep(NA,12),
                            Phylum = rep("Elusimicrobia",12))

Nitrospirae <- data.frame(Taxa = rep("Dummy",12), 
                          Env = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"),
                          Correlation = rep(0,12),
                          Pvalue = rep(0,12),
                          Type = rep("Open woodland",12),
                          AdjPvalue = rep(0,12),
                          Significance = rep(NA,12),
                          Phylum = rep("Nitrospirae",12))

Tenericutes <- data.frame(Taxa = rep("Dummy",12), 
                          Env = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"),
                          Correlation = rep(0,12),
                          Pvalue = rep(0,12),
                          Type = rep("Open woodland",12),
                          AdjPvalue = rep(0,12),
                          Significance = rep(NA,12),
                          Phylum = rep("Tenericutes",12))

WPS <- data.frame(Taxa = rep("Dummy",12), 
                  Env = c("Defol","Temp", "Prec", "pH","Org_C", "Av_P", "NH4_N", "NO3_N", "Ramm", "Rnit", "RH", "Tree_inf"),
                  Correlation = rep(0,12),
                  Pvalue = rep(0,12),
                  Type = rep("Open woodland",12),
                  AdjPvalue = rep(0,12),
                  Significance = rep(NA,12),
                  Phylum = rep("WPS-2",12))

correlacion_nut_func <- rbind(correlacion_nut_func, AD3, Armatimonadetes, Cyanobacteria, Elusimicrobia, Nitrospirae, Tenericutes, WPS)


# Sort phyla according to relative abundance
library(forcats)
## Reordering correlacion_nut_func$Phylum
correlacion_nut_func$Phylum <- factor(correlacion_nut_func$Phylum, levels=c("Acidobacteria", "Actinobacteria", "AD3", "Armatimonadetes", "Bacteroidetes", "Chloroflexi", "Crenarchaeota", "Cyanobacteria", "Elusimicrobia", "Firmicutes", "Gemmatimonadetes", "Nitrospirae", "Planctomycetes", "Proteobacteria", "Tenericutes", "Verrucomicrobia", "WPS-2"))
correlacion_nut_func$Phylum <- fct_rev(correlacion_nut_func$Phylum)

correlacion_nut_func$Type <- "opw"

library(ggplot2)
ggplot(aes(x=Type, y=Phylum, fill=Correlation), data=correlacion_nut_func)+
  geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C", limits=c(-1, 1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)+
  facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x")+
  theme(strip.background = element_rect(fill = "white")) + theme(plot.subtitle = element_text(vjust = 1), 
                                                                 plot.caption = element_text(vjust = 1),
                                                                 text=element_text(family="Times New Roman", face="plain", size=19))
