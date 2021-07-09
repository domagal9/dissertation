#### Replication Materials for Rachel Domagalski's Doctoral Thesis ####
#### Written by Rachel Domagalski, Zachary Neal, Bruce Sagan ####
#### Michigan State University ####
#### Chapter 7 ####

#### Install and Load packages ####
#install.packages("backbone")
#install.packages("igraph")
#install.packages("stringr")
set.seed(19)
library(igraph)  #Used for visualization
library(stringr) #Used for visualization
library(backbone)
sessionInfo()[["otherPkgs"]][["backbone"]][["Version"]]  #Check backbone version, should be 1.5.0

##########################################################
#### Functions to create and visualize igraph objects ####
##########################################################

igraphit <- function(B) {
  #Extract just the positive edges
  B[B<0] <- 0
  
  #Convert it into an igraph object
  iB <- graph_from_adjacency_matrix(B, mode="undirected", diag=FALSE)  #Create igraph object
  
  #Add each senator's political party as a node attribute
  V(iB)$party <- str_sub(gsub(".*\\[(.+)-.*", "\\1",rownames(B)),-2,-2)
  V(iB)$party[V(iB)$party == "I"] <- "D"
  
  #Create a numeric vector for political party
  V(iB)$partynum <- 1  #Republican = 1
  V(iB)$partynum[V(iB)$party == "D"] <- 2 #Democrat = 2
  
  #Define node and edge colors based on political parties
  V(iB)$color[V(iB)$party == "R"] <- rgb(1,0,0,1)
  V(iB)$color[V(iB)$party == "D"] <- rgb(0,0,1,1)
  E(iB)$color <- rgb(.5,0,.5,1)
  republicans <- V(iB)[V(iB)$party=="R"]
  republican.edges <- E(iB)[republicans%--%republicans]
  E(iB)[republican.edges]$color <- rgb(1,0,0,1)
  democrats <- V(iB)[V(iB)$party=="D"]
  democrat.edges <- E(iB)[democrats%--%democrats]
  E(iB)[democrat.edges]$color <- rgb(0,0,1,1)
  
  #Return the compiled igraph object, ready for visualizing
  return(iB)
}

plotit <- function(G) {
  plot(G, 
       vertex.color = V(G)$color,
       vertex.size = 3,
       vertex.frame.color = NA,
       vertex.label.cex = .05,
       edge.color = E(G)$color,
       layout = layout_with_fr)
}

###############################
#### Load example datasets ####
###############################

senate <- read.csv("S114.csv", row.names = 1, header = TRUE)
senate <- as.matrix(senate)
cities <- read.csv(file="https://www.lboro.ac.uk/gawc/datasets/da11.csv", header = TRUE, row.names = 1)
cities <- as.matrix(cities)

###############################
#### Recode cities dataset ####
###############################

cities[cities <= 2] <- 0  #Recode 0s, 1s, and 2s (i.e., Typical and smaller offices) as 0
cities[cities >= 3] <- 1  #Recode 3s, 4s, and 5s (i.e., Larger than typical offices) as 1
cities <- cities[rowSums(cities) != 0,]  #Remove rows that are all 0s (cities with no offices)

#################################
#### Examine Senate Data Set ####
#################################

dim(senate)
senate[1:5, 1:5]
G <- senate%*%t(senate)
dim(G)
G[1:5, 1:2]
G["Booker, C. (NJ-D)", "Warren, E. (MA-D)"]
G["Cruz, T. (TX-R)", "Sanders, B. (VT-I)"]

#################################
#### Examine Cities Data Set ####
#################################

cities[114:117,8:11]  #Look at a portion of the data
rowSums(cities)["AMSTERDAM"]  #Number of firms located in Amsterdam
rowSums(cities)["NEW YORK"]  #Number of firms located in New York
colSums(cities)["KPMG"]  #Number of cities in which KPMG is located
colSums(cities)["HSBC"]  #Number of cities in which HSBC is located
P <- cities %*% t(cities)  #City-to-City projection
P["AMSTERDAM","NEW YORK"]  #Number of firms located in Amsterdam and New York
sort(rowSums(P), decreasing = TRUE)[1:5]  #Report highest degree cities (global network connectivity, GNC)
cor(rowSums(P), rowSums(cities))  #Correlation: Degree and number of firms

##############################
#### Universal Thresholds ####
##############################

### Universal Threshold Set at 0 for the Senate ###
## Extract
universal_bb <- universal(senate, upper = 0, bipartite = TRUE)
universal_bb$backbone[1:5, 1:2]
## Visualize
universal_bb <- universal_bb$backbone
universal_bb.igraph <- igraphit(universal_bb)
plotit(universal_bb.igraph)
## Compute Modularity
universal_bb.mod <- modularity(universal_bb.igraph, V(universal_bb.igraph)$partynum)

### Universal Threshold Set at mean+sd for the Senate ###
## Extract
universal_bb2 <- universal(senate, upper = function(x) mean(x)+sd(x), lower = function(x) mean(x)-sd(x), bipartite = TRUE)
universal_bb2 <- universal_bb2$backbone
## Visualize
universal_bb2.igraph <- igraphit(universal_bb2)
plotit(universal_bb2.igraph)
##Compute Modularity
universal_bb2.mod <- modularity(universal_bb2.igraph, V(universal_bb2.igraph)$partynum)

### Universal Threshold Set at 0 for the Cities ###
## Extract
universal0 <- universal(cities, upper = 0, bipartite = TRUE)  #Compute universal threshold, keeping all non-zero edges
table(universal0$backbone)
## Report Density
mean(universal0$backbone)
sort(rowSums(universal0$backbone), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(universal0$backbone), rowSums(cities))  #Correlation: Degree and number of firms

### Universal Threshold Set at 25 for the Cities ###
## Extract
universal25 <- universal(cities, upper = 25, bipartite = TRUE)  #Compute universal threshold, keeping all non-zero edges
## Report Density
mean(universal25$backbone)
sort(rowSums(universal25$backbone), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(universal25$backbone), rowSums(cities))  #Correlation: Degree and number of firms

### Universal Threshold Set at mean+2sd for the Cities ###
## Extract
universal.meansd <- universal(cities, upper = function(x)mean(x)+2*sd(x), bipartite = TRUE)
## Report Density
mean(universal.meansd$backbone)
sort(rowSums(universal.meansd$backbone), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(universal.meansd$backbone), rowSums(cities))  #Correlation: Degree and number of firms

##########################
#### Fixed Fill Model ####
##########################

### FFM returns NaNs on the Senate Data Set in v1.5.0 ###
# fixedprobs <- fixedfill(senate)

### FFM for Cities ###
fixedprobs <- fixedfill(cities)
fixedbb <- backbone.extract(fixedprobs)
mean(fixedbb)  #Report density
sort(rowSums(fixedbb), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(fixedbb), rowSums(cities))  #Correlation: Degree and number of firms

#########################
#### Fixed Row Model ####
#########################

### FRM for Senate ###
## Extract
fixedrow_probs <- fixedrow(senate)
fixedrow_bb <- backbone.extract(fixedrow_probs, alpha = .01)
## Visualize
fixedrow_bb.igraph <- igraphit(fixedrow_bb)
plotit(fixedrow_bb.igraph)
##Compute Modularity
fixedrow.mod <- modularity(fixedrow_bb.igraph, V(fixedrow_bb.igraph)$partynum)

### FRM on Cities ###
## Extract
rowprobs2 <- fixedrow(cities)
rowbb2 <- backbone.extract(rowprobs2, alpha = .1)
## Report Density
mean(rowbb2)  
sort(rowSums(rowbb2), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(rowbb2), rowSums(cities))  #Correlation: Degree and number of firms

############################
#### Fixed Column Model ####
############################

### FCM for Senate ###
colprobs <- fixedcol(senate)
colbb <- backbone.extract(colprobs, alpha = .01)
## Visualize
fixedcol_bb.igraph <- igraphit(colbb)
plotit(fixedcol_bb.igraph)
##Compute Modularity
fixedcol.mod <- modularity(fixedcol_bb.igraph, V(fixedcol_bb.igraph)$partynum)

### FCM for Cities ###
## Extract
colprobs2 <- fixedcol(cities)
colbb2 <- backbone.extract(colprobs2, alpha = 0.1, signed = FALSE)
## Report Density
mean(colbb2)
sort(rowSums(colbb2), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(colbb2), rowSums(cities))  #Correlation: Degree and number of firms

#################################################
#### Stochastic Degree Sequence Model (SDSM) ####
#################################################

### SDSM for Senate ###
## Extract and time
sdsm <- sdsm(senate)
sdsm_bb <- backbone.extract(sdsm, alpha = .01, narrative = TRUE)

## Visualize
sdsm_bb.igraph <- igraphit(sdsm_bb)
plotit(sdsm_bb.igraph)

##Compute modularity
sdsm.mod <- modularity(sdsm_bb.igraph, V(sdsm_bb.igraph)$partynum)

### SDSM for Cities ###
## Extract
sdsm2 <- sdsm(cities)
sdsmbb2 <- backbone.extract(sdsm2, alpha = 0.2, signed = FALSE, narrative = TRUE)

## Report Density
mean(sdsmbb2)
sort(rowSums(sdsmbb2), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(sdsmbb2), rowSums(cities))  #Correlation: Degree and number of firms

############################################
#### Fixed Degree Sequence Model (FDSM) ####
############################################

### FDSM for Senate ###
## Extract
fdsm <- fdsm(senate, trials = 1000, dyad = c("Booker, C. (NJ-D)", "Warren, E. (MA-D)"), progress = TRUE)
fdsmbb <- backbone.extract(fdsm, signed = TRUE, alpha = 0.01)
## Visualize
fdsm_bb.igraph <- igraphit(fdsmbb)
plotit(fdsmbb.igraph)
## Visualize Histogram
hist(fdsm$dyad_values, 
     freq = FALSE, 
     xlab = "Expected Number of Co-Sponsorships under FDSM", 
     main = NA,
     col = "gray",
     ylim=c(0,.08))
lines(density(fdsm$dyad_values))
## Compute Modularity
fdsm.mod <- modularity(fdsm_bb.igraph, V(fdsm_bb.igraph)$partynum)

### FDSM for Cities ###
## Extract
set.seed(19)
fdsm2 <- fdsm(cities, trials = 10000, progress = TRUE)
fdsmbb2 <- backbone.extract(fdsm2, alpha = 0.1, signed = FALSE)
## Report Density
mean(fdsmbb2)  #Report density
sort(rowSums(fdsmbb2), decreasing = TRUE)[1:5]  #Report highest degree cities
cor(rowSums(fdsmbb2), rowSums(cities))  #Correlation: Degree and number of firms

##################################################
#### Compare FDSM and SDSM in Cities Data Set ####
##################################################

cor(as.vector(fdsmbb2),as.vector(sdsmbb2))  #Correlation: FDSM and SDSM

