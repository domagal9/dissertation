#### Replication Materials for Rachel Domagalski's Doctoral Thesis ####
#### Written by Rachel Domagalski, Zachary Neal, Bruce Sagan ####
#### Michigan State University ####
#### Chapter 8: Study 4 ####

rm(list=ls())
library(backbone)
library(igraph)
set.seed(5)

#### Parameters of bipartite network ####
agents <- 200
artifacts <- 1000
density <- .1
blocks <- 2           #There are two communities
rowdist <- c(1,10)  #Agent degree distributions are right-tailed
coldist <- c(1,10)  #Artifact degree distributions are right-tailed

for (i in 1:10) {  #Loop for replications
  for (block.density in c(seq(from=.5, to=.8, by=.05))) {  #Loop for block density
      message(paste0("Loop = ", i, "   Block Density = ", block.density))

      #Create bipartite, then add block structure
      B <- bipartite.from.distribution(R = agents, C = artifacts, P = density, rowdist = rowdist, coldist = coldist)
      B <- bipartite.add.blocks(B, blocks = blocks, density = block.density)
      
      #Extract backbones
      ffm.bb <- backbone.extract(fixedfill(B), signed = FALSE, alpha = 0.05)
      frm.bb <- backbone.extract(fixedrow(B), signed = FALSE, alpha = 0.05)
      fcm.bb <- backbone.extract(fixedcol(B), signed = FALSE, alpha = 0.05)
      sdsm05.bb <- backbone.extract(sdsm(B), signed = FALSE, alpha = 0.05)
      sdsm13.bb <- backbone.extract(sdsm(B), signed = FALSE, alpha = 0.13)
      fdsm.bb <- backbone.extract(fdsm(B), signed = FALSE, alpha = 0.05)
      
      #Compute modularity
      ffm.i <- graph_from_adjacency_matrix(ffm.bb, mode = "undirected")
      ffm.m <- modularity(ffm.i, as.numeric(as.factor(substr(V(ffm.i)$name,1,1))))
      
      frm.i <- graph_from_adjacency_matrix(frm.bb, mode = "undirected")
      frm.m <- modularity(frm.i, as.numeric(as.factor(substr(V(frm.i)$name,1,1))))
      
      fcm.i <- graph_from_adjacency_matrix(fcm.bb, mode = "undirected")
      fcm.m <- modularity(fcm.i, as.numeric(as.factor(substr(V(fcm.i)$name,1,1))))
      
      sdsm05.i <- graph_from_adjacency_matrix(sdsm05.bb, mode = "undirected")
      sdsm05.m <- modularity(sdsm05.i, as.numeric(as.factor(substr(V(sdsm05.i)$name,1,1))))
      
      sdsm13.i <- graph_from_adjacency_matrix(sdsm13.bb, mode = "undirected")
      sdsm13.m <- modularity(sdsm13.i, as.numeric(as.factor(substr(V(sdsm13.i)$name,1,1))))
      
      fdsm.i <- graph_from_adjacency_matrix(fdsm.bb, mode = "undirected")
      fdsm.m <- modularity(fdsm.i, as.numeric(as.factor(substr(V(fdsm.i)$name,1,1))))
      
      #Post results
      if (exists("modular") == FALSE) {
        modular <- data.frame(loop = i,
                              agents = agents,
                              artifacts = artifacts,
                              density = density,
                              groups = blocks,
                              distrib = paste0(c(rowdist,coldist),collapse=""),
                              block.dens = block.density,
                              ffm.m = ffm.m,
                              frm.m = frm.m,
                              fcm.m = fcm.m,
                              sdsm05.m = sdsm05.m,
                              sdsm13.m = sdsm13.m,
                              fdsm.m = fdsm.m)
      } else modular <- rbind(modular, c(i, agents, artifacts, density, blocks, paste0(c(rowdist,coldist),collapse=""), block.density, ffm.m, frm.m, fcm.m, sdsm05.m, sdsm13.m, fdsm.m))
    }
  }
write.csv(modular, "study4_modularity.csv")

#### Plot ####
rm(list=ls())
library(ggplot2)
library(reshape2)
library(cowplot)
library(latex2exp)
library(backbone)
library(igraph)
library(magick)
set.seed(5)

## Example bipartite
agents <- 200
artifacts <- 1000
density <- .1
blocks <- 2           #There are two communities
rowdist <- c(1,10)  #Agent degree distributions are right-tailed
coldist <- c(1,10)  #Artifact degree distributions are right-tailed

block5 <- bipartite.from.distribution(R = agents, C = artifacts, P = density, rowdist = rowdist, coldist = coldist)
block5 <- bipartite.add.blocks(block5, blocks = blocks, density = .5)
block5 <- block5[order(substr(rownames(block5),1,1),rowSums(block5), decreasing=TRUE), order(substr(colnames(block5),1,1),colSums(block5), decreasing=TRUE)] #Sort by block and degree
block5 <- melt(block5)
block5$value <- factor(block5$value)
block5.plot <-
ggplot(block5, aes(x = Var2, y = Var1, fill = value)) + geom_tile() + #coord_fixed() +
  scale_fill_manual(values = c("white", "black")) + ggtitle("Fraction of within-group edges, W = 0.5") +
  xlab("Group A      Group B\nARTIFACTS") + 
  scale_y_discrete(name = "AGENTS\nGroup B   Group A", limits=rev) +
  theme(axis.title.x=element_text(size = 8), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 8), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 10), legend.position = "none")

block8 <- bipartite.from.distribution(R = agents, C = artifacts, P = density, rowdist = rowdist, coldist = coldist)
block8 <- bipartite.add.blocks(block8, blocks = blocks, density = .8)
block8 <- block8[order(substr(rownames(block8),1,1),rowSums(block8), decreasing=TRUE), order(substr(colnames(block8),1,1),colSums(block8), decreasing=TRUE)] #Sort by block and degree
block8 <- melt(block8)
block8$value <- factor(block8$value)
block8.plot <-
ggplot(block8, aes(x = Var2, y = Var1, fill = value)) + geom_tile() + #coord_fixed() +
  scale_fill_manual(values = c("white", "black")) + ggtitle("Fraction of within-group edges, W = 0.8") +
  xlab("Group A      Group B\nARTIFACTS") + 
  scale_y_discrete(name = "AGENTS\nGroup B   Group A", limits=rev) +
  theme(axis.title.x=element_text(size = 8), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size = 8), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 10), legend.position = "none")

examples <- plot_grid(NULL, block5.plot, NULL, block8.plot, NULL, ncol=1, rel_heights = c(.2, 1, .2, 1, .2))

## Modularity
modular <- read.csv("study4_modularity.csv")
modular <- modular[,c("block.dens", "ffm.m", "frm.m", "fcm.m", "sdsm05.m", "sdsm13.m", "fdsm.m")]
modular <- melt(modular, id.vars = c("block.dens"))
colnames(modular) <- c("density", "model", "modularity")

groups <- 
  ggplot(modular, aes(x=density, y=modularity, color=model)) +
  geom_point(shape = 1) + 
  geom_smooth(method = loess, se = FALSE) +
  scale_x_continuous(name = "Fraction of within-group edges (W)", breaks = c(.5, .55, .6, .65, .7, .75, .8)) +
  ylab(expression("Modularity (Q)")) +
  scale_color_viridis_d(name = "Backbone Model", 
                        limits = c("fdsm.m", "sdsm13.m", "sdsm05.m", "frm.m", "fcm.m", "ffm.m"),
                        labels = c("FDSM", expression(paste("SDSM (",alpha," = 0.13)")), expression(paste("SDSM (",alpha," = 0.05)")), "FRM", "FCM", "FFM")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12), axis.title=element_text(size=12), 
        legend.position = c(0.17, 0.77), legend.text=element_text(size=10), legend.title=element_text(size=12),
        legend.key=element_blank(), legend.text.align = 0)

## Plot backbones
example <- bipartite.from.distribution(R = agents, C = artifacts, P = density, rowdist = rowdist, coldist = coldist)
example <- bipartite.add.blocks(example, blocks = blocks, density = .65)

example.fdsm <- backbone.extract(fdsm(example), signed = FALSE, class = "igraph")
example.fdsm <- delete.vertices(example.fdsm, which(degree(example.fdsm)==0))
V(example.fdsm)$red <- 0
V(example.fdsm)$red[which(substr(V(example.fdsm)$name,1,1)=="A")] <- 1
V(example.fdsm)$blue <- 0
V(example.fdsm)$blue[which(substr(V(example.fdsm)$name,1,1)=="B")] <- 1
png(file="fdsm.png")
plot(example.fdsm,
     vertex.label = NA, vertex.size = 5, edge.color = rgb(0,0,0,.2),
     vertex.color = rgb(V(example.fdsm)$red,0,V(example.fdsm)$blue,.5), vertex.frame.color = NA)
title("FDSM Backbone when W = 0.65",cex.main=2.5)
dev.off()

example.fcm <- backbone.extract(fixedcol(example), signed = FALSE, class = "igraph")
example.fcm <- delete.vertices(example.fcm, which(degree(example.fcm)==0))
V(example.fcm)$red <- 0
V(example.fcm)$red[which(substr(V(example.fcm)$name,1,1)=="A")] <- 1
V(example.fcm)$blue <- 0
V(example.fcm)$blue[which(substr(V(example.fcm)$name,1,1)=="B")] <- 1
png(file="fcm.png")
plot(example.fcm,
     vertex.label = NA, vertex.size = 5, edge.color = rgb(0,0,0,.2),
     vertex.color = rgb(V(example.fcm)$red,0,V(example.fcm)$blue,.5), vertex.frame.color = NA)
title("FCM Backbone when W = 0.65",cex.main=2.5)
dev.off()

fdsm.plot <- ggdraw() + draw_image("fdsm.png")
fcm.plot <- ggdraw() + draw_image("fcm.png")
backbones <- plot_grid(NULL, fdsm.plot, NULL, fcm.plot, NULL, ncol=1, rel_heights = c(.2, 1, .2, 1, .2))
unlink("fdsm.png")
unlink("fcm.png")

## Combine
plot_grid(examples, NULL, groups, NULL, backbones, ncol=5, labels=c("A","","B","","C"), rel_widths = c(1, .1, 2, .1, 1))
ggsave("study4_plot.pdf", width = 12, height = 5, device = "pdf")
