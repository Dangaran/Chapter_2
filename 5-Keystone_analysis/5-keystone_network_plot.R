bact_network_deh_igraph
keystones_bact_for_mod$ID

# Shapes
shapes <- setdiff(shapes(), "")
shapes

mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

shapes("circle")$clip <- shapes("circle")$clip

# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}
# no clipping, edges will be below the vertices anyway
add_shape("star", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=5))


# Bacteria Forest ----
# Loop to find number of keystones
keyst_for_num <-c()
for(i in 1:nrow(keystones_bact_for_mod)){
  keyst_for_num[i] <- print(which(V(bact_network_for_igraph)$name == keystones_bact_for_mod$ID[i]))
}

# List with circundant nodes
library(igraph)
combine_for_keyst_list <- make_ego_graph(bact_network_for_igraph, order = 1, nodes = keyst_for_num, 
               mode = c("all"), mindist = 0)

# Combine list in a single graph https://igraph.org/r/doc/union.igraph.html
combine_for_keyst <- print_all(combine_for_keyst_list[[1]] %u% combine_for_keyst_list[[2]]%u% 
                 combine_for_keyst_list[[3]] %u% combine_for_keyst_list[[4]] %u%
                 combine_for_keyst_list[[5]] %u% combine_for_keyst_list[[6]] %u%
                 combine_for_keyst_list[[7]] %u% combine_for_keyst_list[[8]] %u%
                 combine_for_keyst_list[[9]] %u% combine_for_keyst_list[[10]] %u%
                 combine_for_keyst_list[[11]] %u% combine_for_keyst_list[[12]] %u%
                 combine_for_keyst_list[[13]] %u% combine_for_keyst_list[[14]] %u%
                 combine_for_keyst_list[[15]] %u% combine_for_keyst_list[[16]] %u%
                 combine_for_keyst_list[[17]] %u% combine_for_keyst_list[[18]] %u%
                 combine_for_keyst_list[[19]] %u% combine_for_keyst_list[[20]] %u%
                 combine_for_keyst_list[[21]])

V(combine_for_keyst)$shape <- "circle"

# Change color to white
V(combine_for_keyst)$color <- "white"

# Change size
V(combine_for_keyst)$size <- 1


# Extract keystone type
Deg_for <- subset(keystones_bact_for_mod, Tipo_Keyst == "Degree")
BTW_for <- subset(keystones_bact_for_mod, Tipo_Keyst == "BTW")
Deg_BTW_for <- subset(keystones_bact_for_mod, Tipo_Keyst == "Degree_BTW")
unique(keystones_bact_for_mod$Phylum)
Acidobacteria_for <- subset(keystones_bact_for_mod, Phylum == "Acidobacteria")
Actinobacteria_for <- subset(keystones_bact_for_mod, Phylum == "Actinobacteria")
Bacteroidetes_for <- subset(keystones_bact_for_mod, Phylum == "Bacteroidetes")
Chloroflexi_for <- subset(keystones_bact_for_mod, Phylum == "Chloroflexi")
Proteobacteria_for <- subset(keystones_bact_for_mod, Phylum == "Proteobacteria")
Verrucomicrobia_for <- subset(keystones_bact_for_mod, Phylum == "Verrucomicrobia")


# Modify shape and color for keystones
for (i in 1:nrow(keystones_bact_for_mod)) {
  V(combine_for_keyst)$size[which(V(combine_for_keyst)$name == Deg_for$ID[i])] <- 5
  V(combine_for_keyst)$size[which(V(combine_for_keyst)$name == BTW_for$ID[i])] <- 7
  V(combine_for_keyst)$size[which(V(combine_for_keyst)$name == Deg_BTW_for$ID[i])] <- 7
  V(combine_for_keyst)$shape[which(V(combine_for_keyst)$name == Deg_for$ID[i])] <- "square"
  V(combine_for_keyst)$shape[which(V(combine_for_keyst)$name == BTW_for$ID[i])] <- "triangle"
  V(combine_for_keyst)$shape[which(V(combine_for_keyst)$name == Deg_BTW_for$ID[i])] <- "star"
  V(combine_for_keyst)$color[which(V(combine_for_keyst)$name == Acidobacteria_for$ID[i])] <- "#0000CC"
  V(combine_for_keyst)$color[which(V(combine_for_keyst)$name == Actinobacteria_for$ID[i])] <- "#666666"
  V(combine_for_keyst)$color[which(V(combine_for_keyst)$name == Bacteroidetes_for$ID[i])] <- "#339999"
  V(combine_for_keyst)$color[which(V(combine_for_keyst)$name == Chloroflexi_for$ID[i])] <- "#336600"
  V(combine_for_keyst)$color[which(V(combine_for_keyst)$name == Proteobacteria_for$ID[i])] <- "orange"
  V(combine_for_keyst)$color[which(V(combine_for_keyst)$name == Verrucomicrobia_for$ID[i])] <- "#CC0033"
  
}

# Plotting
plot(combine_for_keyst, layout=layout_with_graphopt, vertex.size=combine_for_keyst$size, vertex.label=NA)




# Bacteria Dehesa ----
# Loop to find number of keystones
keyst_deh_num <-c()
for(i in 1:nrow(keystones_bact_deh_mod)){
  keyst_deh_num[i] <- print(which(V(bact_network_deh_igraph)$name == keystones_bact_deh_mod$ID[i]))
}

# List with circundant nodes
library(igraph)
combine_deh_keyst_list <- make_ego_graph(bact_network_deh_igraph, order = 1, nodes = keyst_deh_num, 
                                         mode = c("all"), mindist = 0)

# Combine list in a single graph https://igraph.org/r/doc/union.igraph.html
combine_deh_keyst <- print_all(combine_deh_keyst_list[[1]] %u% combine_deh_keyst_list[[2]]%u% 
                                 combine_deh_keyst_list[[3]] %u% combine_deh_keyst_list[[4]] %u%
                                 combine_deh_keyst_list[[5]] %u% combine_deh_keyst_list[[6]] %u%
                                 combine_deh_keyst_list[[7]] %u% combine_deh_keyst_list[[8]] %u%
                                 combine_deh_keyst_list[[9]] %u% combine_deh_keyst_list[[10]] %u%
                                 combine_deh_keyst_list[[11]] %u% combine_deh_keyst_list[[12]] %u%
                                 combine_deh_keyst_list[[13]])

V(combine_deh_keyst)$shape <- "circle"

# Change color to white
V(combine_deh_keyst)$color <- "white"

# Change size
V(combine_deh_keyst)$size <- 1

# Extract keystone type
Deg_deh <- subset(keystones_bact_deh_mod, Tipo_Keyst == "Degree")
BTW_deh <- subset(keystones_bact_deh_mod, Tipo_Keyst == "BTW")
Deg_BTW_deh <- subset(keystones_bact_deh_mod, Tipo_Keyst == "Degree_BTW")
unique(keystones_bact_deh_mod$Phylum)
Acidobacteria_deh <- subset(keystones_bact_deh_mod, Phylum == "Acidobacteria")
Bacteroidetes_deh <- subset(keystones_bact_deh_mod, Phylum == "Bacteroidetes")
Chloroflexi_deh <- subset(keystones_bact_deh_mod, Phylum == "Chloroflexi")
Proteobacteria_deh <- subset(keystones_bact_deh_mod, Phylum == "Proteobacteria")
Gemmatimonadetes_deh <- subset(keystones_bact_deh_mod, Phylum == "Gemmatimonadetes")


# Modify shape and color for keystones
for (i in 1:nrow(keystones_bact_deh_mod)) {
  V(combine_deh_keyst)$size[which(V(combine_deh_keyst)$name == Deg_deh$ID[i])] <- 5
  V(combine_deh_keyst)$size[which(V(combine_deh_keyst)$name == BTW_deh$ID[i])] <- 7
  V(combine_deh_keyst)$size[which(V(combine_deh_keyst)$name == Deg_BTW_deh$ID[i])] <- 7
  V(combine_deh_keyst)$shape[which(V(combine_deh_keyst)$name == Deg_deh$ID[i])] <- "square"
  V(combine_deh_keyst)$shape[which(V(combine_deh_keyst)$name == BTW_deh$ID[i])] <- "triangle"
  V(combine_deh_keyst)$shape[which(V(combine_deh_keyst)$name == Deg_BTW_deh$ID[i])] <- "star"
  V(combine_deh_keyst)$color[which(V(combine_deh_keyst)$name == Acidobacteria_deh$ID[i])] <- "#0000CC"
  V(combine_deh_keyst)$color[which(V(combine_deh_keyst)$name == Bacteroidetes_deh$ID[i])] <- "#339999"
  V(combine_deh_keyst)$color[which(V(combine_deh_keyst)$name == Chloroflexi_deh$ID[i])] <- "#336600"
  V(combine_deh_keyst)$color[which(V(combine_deh_keyst)$name == Proteobacteria_deh$ID[i])] <- "orange"
  V(combine_deh_keyst)$color[which(V(combine_deh_keyst)$name == Gemmatimonadetes_deh$ID[i])] <- "#990033"
  
}


# Plotting
plot(combine_deh_keyst, layout=layout_with_graphopt, vertex.size=combine_deh_keyst$size, vertex.label=NA)



# Bacteria Open woodland ----
# Loop to find number of keystones
keyst_opw_num <-c()
for(i in 1:nrow(keystones_bact_opw_mod)){
  keyst_opw_num[i] <- print(which(V(bact_network_opw_igraph)$name == keystones_bact_opw_mod$ID[i]))
}

# List with circundant nodes
library(igraph)
combine_opw_keyst_list <- make_ego_graph(bact_network_opw_igraph, order = 1, nodes = keyst_opw_num, 
                                         mode = c("all"), mindist = 0)

# Combine list in a single graph https://igraph.org/r/doc/union.igraph.html
combine_opw_keyst <- print_all(combine_opw_keyst_list[[1]] %u% combine_opw_keyst_list[[2]]%u% 
                                 combine_opw_keyst_list[[3]] %u% combine_opw_keyst_list[[4]] %u%
                                 combine_opw_keyst_list[[5]] %u% combine_opw_keyst_list[[6]] %u%
                                 combine_opw_keyst_list[[7]] %u% combine_opw_keyst_list[[8]] %u%
                                 combine_opw_keyst_list[[9]] %u% combine_opw_keyst_list[[10]] %u%
                                 combine_opw_keyst_list[[11]] %u% combine_opw_keyst_list[[12]] %u%
                                 combine_opw_keyst_list[[13]] %u% combine_opw_keyst_list[[14]] %u%
                                 combine_opw_keyst_list[[15]] %u% combine_opw_keyst_list[[16]] %u%
                                 combine_opw_keyst_list[[17]] %u% combine_opw_keyst_list[[18]] %u%
                                 combine_opw_keyst_list[[19]])

V(combine_opw_keyst)$shape <- "circle"

# Change color to white
V(combine_opw_keyst)$color <- "white"

# Change size
V(combine_opw_keyst)$size <- 1

# Extract keystone type
Deg_opw <- subset(keystones_bact_opw_mod, Tipo_Keyst == "Degree")
BTW_opw <- subset(keystones_bact_opw_mod, Tipo_Keyst == "BTW")
Deg_BTW_opw <- subset(keystones_bact_opw_mod, Tipo_Keyst == "Degree_BTW")
unique(keystones_bact_opw_mod$Phylum)
Acidobacteria_opw <- subset(keystones_bact_opw_mod, Phylum == "Acidobacteria")
Actinobacteria_opw <- subset(keystones_bact_opw_mod, Phylum == "Actinobacteria")
Bacteroidetes_opw <- subset(keystones_bact_opw_mod, Phylum == "Bacteroidetes")
Chloroflexi_opw <- subset(keystones_bact_opw_mod, Phylum == "Chloroflexi")
Proteobacteria_opw <- subset(keystones_bact_opw_mod, Phylum == "Proteobacteria")
Verrucomicrobia_opw <- subset(keystones_bact_opw_mod, Phylum == "Verrucomicrobia")
Planctomycetes_opw <- subset(keystones_bact_opw_mod, Phylum == "Planctomycetes")


# Modify shape and color for keystones
for (i in 1:nrow(keystones_bact_opw_mod)) {
  V(combine_opw_keyst)$size[which(V(combine_opw_keyst)$name == Deg_opw$ID[i])] <- 5
  V(combine_opw_keyst)$size[which(V(combine_opw_keyst)$name == BTW_opw$ID[i])] <- 7
  V(combine_opw_keyst)$size[which(V(combine_opw_keyst)$name == Deg_BTW_opw$ID[i])] <- 7
  V(combine_opw_keyst)$shape[which(V(combine_opw_keyst)$name == Deg_opw$ID[i])] <- "square"
  V(combine_opw_keyst)$shape[which(V(combine_opw_keyst)$name == BTW_opw$ID[i])] <- "triangle"
  V(combine_opw_keyst)$shape[which(V(combine_opw_keyst)$name == Deg_BTW_opw$ID[i])] <- "star"
  V(combine_opw_keyst)$color[which(V(combine_opw_keyst)$name == Acidobacteria_opw$ID[i])] <- "#0000CC"
  V(combine_opw_keyst)$color[which(V(combine_opw_keyst)$name == Actinobacteria_opw$ID[i])] <- "#666666"
  V(combine_opw_keyst)$color[which(V(combine_opw_keyst)$name == Bacteroidetes_opw$ID[i])] <- "#339999"
  V(combine_opw_keyst)$color[which(V(combine_opw_keyst)$name == Chloroflexi_opw$ID[i])] <- "#336600"
  V(combine_opw_keyst)$color[which(V(combine_opw_keyst)$name == Proteobacteria_opw$ID[i])] <- "orange"
  V(combine_opw_keyst)$color[which(V(combine_opw_keyst)$name == Verrucomicrobia_opw$ID[i])] <- "#CC0033"
  V(combine_opw_keyst)$color[which(V(combine_opw_keyst)$name == Planctomycetes_opw$ID[i])] <- "yellow"
  
}


# Plotting
plot(combine_opw_keyst, layout=layout_with_graphopt, vertex.size=combine_opw_keyst$size, vertex.label=NA)




# Fungi Forests ----
# Loop to find number of keystones
keyst_fungi_for_num <-c()
for(i in 1:nrow(keystones_fungi_for_mod)){
  keyst_fungi_for_num[i] <- print(which(V(fungi_network_for_igraph)$name == keystones_fungi_for_mod$ID[i]))
}

# List with circundant nodes
library(igraph)
combine_fungi_for_keyst_list <- make_ego_graph(fungi_network_for_igraph, order = 1, nodes = keyst_fungi_for_num, 
                                         mode = c("all"), mindist = 0)

# Combine list in a single graph https://igraph.org/r/doc/union.igraph.html
combine_fungi_for_keyst <- print_all(combine_fungi_for_keyst_list[[1]] %u% combine_fungi_for_keyst_list[[2]]%u% 
                                 combine_fungi_for_keyst_list[[3]] %u% combine_fungi_for_keyst_list[[4]] %u%
                                 combine_fungi_for_keyst_list[[5]] %u% combine_fungi_for_keyst_list[[6]] %u%
                                 combine_fungi_for_keyst_list[[7]] %u% combine_fungi_for_keyst_list[[8]] %u%
                                 combine_fungi_for_keyst_list[[9]] %u% combine_fungi_for_keyst_list[[10]] %u%
                                 combine_fungi_for_keyst_list[[11]] %u% combine_fungi_for_keyst_list[[12]] %u%
                                 combine_fungi_for_keyst_list[[13]] %u% combine_fungi_for_keyst_list[[14]] %u%
                                 combine_fungi_for_keyst_list[[15]] %u% combine_fungi_for_keyst_list[[16]] %u%
                                 combine_fungi_for_keyst_list[[17]] %u% combine_fungi_for_keyst_list[[18]] %u%
                                 combine_fungi_for_keyst_list[[19]])

V(combine_fungi_for_keyst)$shape <- "circle"

# Change color to white
V(combine_fungi_for_keyst)$color <- "white"

# Change size
V(combine_fungi_for_keyst)$size <- 1

# Extract keystone type
Deg_fungi_for <- subset(keystones_fungi_for_mod, Tipo_Keyst == "Degree")
BTW_fungi_for <- subset(keystones_fungi_for_mod, Tipo_Keyst == "BTW")
Deg_BTW_fungi_for <- subset(keystones_fungi_for_mod, Tipo_Keyst == "Degree_BTW")
unique(keystones_fungi_for_mod$Phylum)
Ascomycota_for <- subset(keystones_fungi_for_mod, Phylum == "p__Ascomycota")
Basidiomycota_for <- subset(keystones_fungi_for_mod, Phylum == "p__Basidiomycota")


# Modify shape and color for keystones
for (i in 1:nrow(keystones_fungi_for_mod)) {
  V(combine_fungi_for_keyst)$size[which(V(combine_fungi_for_keyst)$name == Deg_fungi_for$ID[i])] <- 5
  V(combine_fungi_for_keyst)$size[which(V(combine_fungi_for_keyst)$name == BTW_fungi_for$ID[i])] <- 7
  V(combine_fungi_for_keyst)$size[which(V(combine_fungi_for_keyst)$name == Deg_BTW_fungi_for$ID[i])] <- 7
  V(combine_fungi_for_keyst)$shape[which(V(combine_fungi_for_keyst)$name == Deg_fungi_for$ID[i])] <- "square"
  V(combine_fungi_for_keyst)$shape[which(V(combine_fungi_for_keyst)$name == BTW_fungi_for$ID[i])] <- "triangle"
  V(combine_fungi_for_keyst)$shape[which(V(combine_fungi_for_keyst)$name == Deg_BTW_fungi_for$ID[i])] <- "star"
  V(combine_fungi_for_keyst)$color[which(V(combine_fungi_for_keyst)$name == Ascomycota_for$ID[i])] <- "#009900"
  V(combine_fungi_for_keyst)$color[which(V(combine_fungi_for_keyst)$name == Basidiomycota_for$ID[i])] <- "orange"
}

# Plotting
plot(combine_fungi_for_keyst, layout=layout_with_graphopt, vertex.size=combine_fungi_for_keyst$size, vertex.label=NA)




# Fungi Dehesas ----
# Loop to find number of keystones
keyst_fungi_deh_num <-c()
for(i in 1:nrow(keystones_fungi_deh_mod)){
  keyst_fungi_deh_num[i] <- print(which(V(fungi_network_deh_igraph)$name == keystones_fungi_deh_mod$ID[i]))
}

# List with circundant nodes
library(igraph)
combine_fungi_deh_keyst_list <- make_ego_graph(fungi_network_deh_igraph, order = 1, nodes = keyst_fungi_deh_num, 
                                               mode = c("all"), mindist = 0)

# Combine list in a single graph https://igraph.org/r/doc/union.igraph.html
combine_fungi_deh_keyst <- print_all(combine_fungi_deh_keyst_list[[1]] %u% combine_fungi_deh_keyst_list[[2]]%u% 
                                       combine_fungi_deh_keyst_list[[3]] %u% combine_fungi_deh_keyst_list[[4]] %u%
                                       combine_fungi_deh_keyst_list[[5]] %u% combine_fungi_deh_keyst_list[[6]] %u%
                                       combine_fungi_deh_keyst_list[[7]] %u% combine_fungi_deh_keyst_list[[8]] %u%
                                       combine_fungi_deh_keyst_list[[9]] %u% combine_fungi_deh_keyst_list[[10]] %u%
                                       combine_fungi_deh_keyst_list[[11]] %u% combine_fungi_deh_keyst_list[[12]] %u%
                                       combine_fungi_deh_keyst_list[[13]] %u% combine_fungi_deh_keyst_list[[14]] %u%
                                       combine_fungi_deh_keyst_list[[15]] %u% combine_fungi_deh_keyst_list[[16]] %u%
                                       combine_fungi_deh_keyst_list[[17]])

V(combine_fungi_deh_keyst)$shape <- "circle"

# Change color to white
V(combine_fungi_deh_keyst)$color <- "white"

# Change size
V(combine_fungi_deh_keyst)$size <- 1

# Extract keystone type
Deg_fungi_deh <- subset(keystones_fungi_deh_mod, Tipo_Keyst == "Degree")
BTW_fungi_deh <- subset(keystones_fungi_deh_mod, Tipo_Keyst == "BTW")
Deg_BTW_fungi_deh <- subset(keystones_fungi_deh_mod, Tipo_Keyst == "Degree_BTW")
unique(keystones_fungi_deh_mod$Phylum)
Ascomycota_deh <- subset(keystones_fungi_deh_mod, Phylum == "p__Ascomycota")
Basidiomycota_deh <- subset(keystones_fungi_deh_mod, Phylum == "p__Basidiomycota")


# Modify shape and color for keystones
for (i in 1:nrow(keystones_fungi_deh_mod)) {
  V(combine_fungi_deh_keyst)$size[which(V(combine_fungi_deh_keyst)$name == Deg_fungi_deh$ID[i])] <- 5
  V(combine_fungi_deh_keyst)$size[which(V(combine_fungi_deh_keyst)$name == BTW_fungi_deh$ID[i])] <- 7
  V(combine_fungi_deh_keyst)$size[which(V(combine_fungi_deh_keyst)$name == Deg_BTW_fungi_deh$ID[i])] <- 7
  V(combine_fungi_deh_keyst)$shape[which(V(combine_fungi_deh_keyst)$name == Deg_fungi_deh$ID[i])] <- "square"
  V(combine_fungi_deh_keyst)$shape[which(V(combine_fungi_deh_keyst)$name == BTW_fungi_deh$ID[i])] <- "triangle"
  V(combine_fungi_deh_keyst)$shape[which(V(combine_fungi_deh_keyst)$name == Deg_BTW_fungi_deh$ID[i])] <- "star"
  V(combine_fungi_deh_keyst)$color[which(V(combine_fungi_deh_keyst)$name == Ascomycota_deh$ID[i])] <- "#009900"
  V(combine_fungi_deh_keyst)$color[which(V(combine_fungi_deh_keyst)$name == Basidiomycota_deh$ID[i])] <- "orange"
}


# Plotting
plot(combine_fungi_deh_keyst, layout=layout_with_graphopt, vertex.size=combine_fungi_deh_keyst$size, vertex.label=NA)


# Fungi Open woodland ----
# Loop to find number of keystones
keyst_fungi_opw_num <-c()
for(i in 1:nrow(keystones_fungi_opw_mod)){
  keyst_fungi_opw_num[i] <- print(which(V(fungi_network_opw_igraph)$name == keystones_fungi_opw_mod$ID[i]))
}

# List with circundant nodes
library(igraph)
combine_fungi_opw_keyst_list <- make_ego_graph(fungi_network_opw_igraph, order = 1, nodes = keyst_fungi_opw_num, 
                                               mode = c("all"), mindist = 0)

# Combine list in a single graph https://igraph.org/r/doc/union.igraph.html
combine_fungi_opw_keyst <- print_all(combine_fungi_opw_keyst_list[[1]] %u% combine_fungi_opw_keyst_list[[2]]%u% 
                                       combine_fungi_opw_keyst_list[[3]] %u% combine_fungi_opw_keyst_list[[4]] %u%
                                       combine_fungi_opw_keyst_list[[5]] %u% combine_fungi_opw_keyst_list[[6]] %u%
                                       combine_fungi_opw_keyst_list[[7]] %u% combine_fungi_opw_keyst_list[[8]] %u%
                                       combine_fungi_opw_keyst_list[[9]] %u% combine_fungi_opw_keyst_list[[10]] %u%
                                       combine_fungi_opw_keyst_list[[11]] %u% combine_fungi_opw_keyst_list[[12]] %u%
                                       combine_fungi_opw_keyst_list[[13]] %u% combine_fungi_opw_keyst_list[[14]] %u%
                                       combine_fungi_opw_keyst_list[[15]] %u% combine_fungi_opw_keyst_list[[16]] %u%
                                       combine_fungi_opw_keyst_list[[17]] %u% combine_fungi_opw_keyst_list[[18]] %u%
                                       combine_fungi_opw_keyst_list[[19]] %u% combine_fungi_opw_keyst_list[[20]] %u%
                                       combine_fungi_opw_keyst_list[[21]] %u% combine_fungi_opw_keyst_list[[22]] %u%
                                       combine_fungi_opw_keyst_list[[23]] %u% combine_fungi_opw_keyst_list[[24]])

V(combine_fungi_opw_keyst)$shape <- "circle"

# Change color to white
V(combine_fungi_opw_keyst)$color <- "white"

# Change size
V(combine_fungi_opw_keyst)$size <- 1

# Extract keystone type
Deg_fungi_opw <- subset(keystones_fungi_opw_mod, Tipo_Keyst == "Degree")
BTW_fungi_opw <- subset(keystones_fungi_opw_mod, Tipo_Keyst == "BTW")
Deg_BTW_fungi_opw <- subset(keystones_fungi_opw_mod, Tipo_Keyst == "Degree_BTW")
unique(keystones_fungi_opw_mod$Phylum)
Ascomycota_opw <- subset(keystones_fungi_opw_mod, Phylum == "p__Ascomycota")
Basidiomycota_opw <- subset(keystones_fungi_opw_mod, Phylum == "p__Basidiomycota")
Mortierellomycota_opw <- subset(keystones_fungi_opw_mod, Phylum == "p__Mortierellomycota")
Mucoromycota_opw <- subset(keystones_fungi_opw_mod, Phylum == "p__Mucoromycota")


# Modify shape and color for keystones
for (i in 1:nrow(keystones_fungi_opw_mod)) {
  V(combine_fungi_opw_keyst)$size[which(V(combine_fungi_opw_keyst)$name == Deg_fungi_opw$ID[i])] <- 5
  V(combine_fungi_opw_keyst)$size[which(V(combine_fungi_opw_keyst)$name == BTW_fungi_opw$ID[i])] <- 7
  V(combine_fungi_opw_keyst)$size[which(V(combine_fungi_opw_keyst)$name == Deg_BTW_fungi_opw$ID[i])] <- 7
  V(combine_fungi_opw_keyst)$shape[which(V(combine_fungi_opw_keyst)$name == Deg_fungi_opw$ID[i])] <- "square"
  V(combine_fungi_opw_keyst)$shape[which(V(combine_fungi_opw_keyst)$name == BTW_fungi_opw$ID[i])] <- "triangle"
  V(combine_fungi_opw_keyst)$shape[which(V(combine_fungi_opw_keyst)$name == Deg_BTW_fungi_opw$ID[i])] <- "star"
  V(combine_fungi_opw_keyst)$color[which(V(combine_fungi_opw_keyst)$name == Ascomycota_opw$ID[i])] <- "#009900"
  V(combine_fungi_opw_keyst)$color[which(V(combine_fungi_opw_keyst)$name == Basidiomycota_opw$ID[i])] <- "orange"
  V(combine_fungi_opw_keyst)$color[which(V(combine_fungi_opw_keyst)$name == Mortierellomycota_opw$ID[i])] <- "#0000CC"
  V(combine_fungi_opw_keyst)$color[which(V(combine_fungi_opw_keyst)$name == Mucoromycota_opw$ID[i])] <- "#990000"
  
}


# Plotting
plot(combine_fungi_opw_keyst, layout=layout_with_graphopt, vertex.size=combine_fungi_opw_keyst$size, vertex.label=NA)
legend(1.2,0.5, legend=dummy_fungi$Phylum,pch=21, pt.bg = dummy_fungi$color, bg=NULL, title="Phylum", box.lty=0)

dummy_fungi$color <- c("red", "orange", "violet", "darkgreen")
dummy_fungi$Phylum <- c("p__Ascomycota", "p__Basidiomycota", "p__Mortierellomycota", "p__Mucoromycota")


par(mfrow=c(2,3))
plot(combine_for_keyst, layout=layout_with_graphopt, vertex.size=combine_for_keyst$size, vertex.label=NA)
plot(combine_deh_keyst, layout=layout_with_graphopt, vertex.size=combine_deh_keyst$size, vertex.label=NA)
plot(combine_opw_keyst, layout=layout_with_graphopt, vertex.size=combine_opw_keyst$size, vertex.label=NA)
plot(combine_fungi_for_keyst, layout=layout_with_graphopt, vertex.size=combine_fungi_for_keyst$size, vertex.label=NA)
plot(combine_fungi_deh_keyst, layout=layout_with_graphopt, vertex.size=combine_fungi_deh_keyst$size, vertex.label=NA)
plot(combine_fungi_opw_keyst, layout=layout_with_graphopt, vertex.size=combine_fungi_opw_keyst$size, vertex.label=NA)
legend(1.2,0.5, legend=dummy_fungi$Phylum,pch=21, pt.bg = dummy_fungi$color, bg=NULL, title="Phylum", box.lty=0)
dev.off()
