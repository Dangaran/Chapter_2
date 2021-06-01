##%######################################################%##
#                                                          #
####               Venn diagram BACTERIA                ####
#                                                          #
##%######################################################%##
keystones_bact_for_mod <- read.csv("Nuevas_keystones/keystones_bact_for_mod_mod.csv", h = T, sep = ";")
keystones_bact_deh_mod <- read.csv("Nuevas_keystones/keystones_bact_deh_mod.csv", h = T, sep = ";")
keystones_bact_opw_mod <- read.csv("Nuevas_keystones/keystones_bact_opw_mod.csv", h = T, sep = ";")

keystones_fungi_for_mod <- read.csv("Nuevas_keystones/keystones_fungi_for_mod.csv", h = T, sep = ";")
keystones_fungi_deh_mod <- read.csv("Nuevas_keystones/keystones_fungi_deh_mod.csv", h = T, sep = ";")
keystones_fungi_opw_mod <- read.csv("Nuevas_keystones/keystones_fungi_opw_mod.csv", h = T, sep = ";")

library(VennDiagram)
nrow(keystones_bact_for_mod)
nrow(keystones_bact_deh_mod)
nrow(keystones_bact_opw_mod)

grid.newpage()
venn.plot_bact <- draw.triple.venn(
  area1 = 21,
  area2 = 19,
  area3 = 13,
  n12 = 5,
  n23 = 0,
  n13 = 3,
  n123 = 0,
  category = c("Forest", "Open woodland", "Dehesa"),
  #fill = c("darkolivegreen3", "yellowgreen", "burlywood3"),
  lty = c(1,2,3),
  cex = 2,
  cat.cex = 2)
grid.draw(venn.plot_bact)
grid.newpage()



# Forest vs Dehesa
data.frame(Function=keystones_bact_deh_mod[match(keystones_bact_for_mod$ID, keystones_bact_deh_mod$ID), 1]) 
grid.newpage()
draw.pairwise.venn(area1 = 21, area2 = 13, cross.area = 3, category = c("Forest", "Dehesa"), lty = "blank", 
                   fill = c("darkolivegreen3", "burlywood3"))

# Forest vs Open woodland
data.frame(Function=keystones_bact_opw_mod[match(keystones_bact_for_mod$ID, keystones_bact_opw_mod$ID), 1]) 
grid.newpage()
draw.pairwise.venn(area1 = 21, area2 = 19, cross.area = 5, category = c("Forest", "Open woodland"), lty = "blank", 
                   fill = c("darkolivegreen3", "burlywood3"))

# Dehesa vs Open woodland
data.frame(Function=keystones_bact_deh_mod[match(keystones_bact_opw_mod$ID, keystones_bact_deh_mod$ID), 1]) 
grid.newpage()
draw.pairwise.venn(area1 = 13, area2 = 19, cross.area = 0, category = c("Dehesa", "Open woodland"), lty = "blank", 
                 fill = c("darkolivegreen3", "burlywood3"))


grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=3, widths = unit(rep(1/3,3), "npc"))))
pushViewport(viewport(layout.pos.row =1))
draw.pairwise.venn(area1 = 21, area2 = 13, cross.area = 3, category = c("Forest", "Dehesa"), lty = "blank", fill = c("darkolivegreen3", "burlywood3"))
popViewport()
pushViewport(viewport(layout.pos.row=2))
draw.pairwise.venn(area1 = 21, area2 = 19, cross.area = 5, category = c("Forest", "Open woodland"), lty = "blank", fill = c("darkolivegreen3", "yellowgreen"))
popViewport()
pushViewport(viewport(layout.pos.row=3))
draw.pairwise.venn(area1 = 13, area2 = 19, cross.area = 0, category = c("Dehesa", "Open woodland"), lty = "blank", fill = c("burlywood3", "yellowgreen"))
popViewport(0)

##%######################################################%##
#                                                          #
####                Venn diagram Fungi                  ####
#                                                          #
##%######################################################%##
library(VennDiagram)
nrow(keystones_fungi_for_mod)
nrow(keystones_fungi_deh_mod)
nrow(keystones_fungi_opw_mod)

grid.newpage()
venn.plot_fungi <- draw.triple.venn(
  area1 = 19,
  area2 = 24,
  area3 = 17,
  n12 = 0,
  n23 = 0,
  n13 = 1,
  n123 = 0,
  category = c("Forest", "Open woodland", "Dehesa"),
  #fill = c("darkolivegreen3", "yellowgreen", "burlywood3"),
  lty = c(1,2,3),
  cex = 2,
  cat.cex = 2)
grid.draw(venn.plot_fungi)
grid.newpage()


# Forest vs Dehesa
data.frame(Function=keystones_fungi_deh_mod[match(keystones_fungi_for_mod$ID, keystones_fungi_deh_mod$ID), 1])
grid.newpage()
draw.pairwise.venn(area1 = 19, area2 = 17, cross.area = 1, category = c("Forest", "Dehesa"), lty = "blank", 
                   fill = c("darkolivegreen3", "burlywood3"))

# Forest vs Open woodland
data.frame(Function=keystones_fungi_for_mod[match(keystones_fungi_opw_mod$ID, keystones_fungi_for_mod$ID), 1])
grid.newpage()
draw.pairwise.venn(area1 = 19, area2 = 24, cross.area = 0, category = c("Forest", "Open woodland"), lty = "blank", 
                   fill = c("darkolivegreen3", "burlywood3"))

# Dehesa vs Open woodland
data.frame(Function=keystones_fungi_deh_mod[match(keystones_fungi_opw_mod$ID, keystones_fungi_deh_mod$ID), 1])
grid.newpage()
draw.pairwise.venn(area1 = 17, area2 = 24, cross.area = 0, category = c("Dehesa", "Open woodland"), lty = "blank", 
                   fill = c("darkolivegreen3", "burlywood3"))


grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=3, widths = unit(rep(1/3,3), "npc"))))
pushViewport(viewport(layout.pos.row =1))
draw.pairwise.venn(area1 = 19, area2 = 17, cross.area = 1, category = c("Forest", "Dehesa"), lty = "blank", fill = c("darkolivegreen3", "burlywood3"))
popViewport()
pushViewport(viewport(layout.pos.row=2))
draw.pairwise.venn(area1 = 19, area2 = 24, cross.area = 0, category = c("Forest", "Open woodland"), lty = "blank", fill = c("darkolivegreen3", "yellowgreen"))
popViewport()
pushViewport(viewport(layout.pos.row=3))
draw.pairwise.venn(area1 = 17, area2 = 24, cross.area = 0, category = c("Dehesa", "Open woodland"), lty = "blank", fill = c("burlywood3", "yellowgreen"))
popViewport(0)
