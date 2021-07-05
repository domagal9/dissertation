#### Replication Materials for Rachel Domagalski's Doctoral Thesis ####
#### Written by Rachel Domagalski, Zachary Neal, Bruce Sagan ####
#### Michigan State University ####
#### Chapter 8: Study 3 ####

#### Load packages and set seed ####
rm(list=ls())
library(backbone)
set.seed(5)

## Define jaccard coefficient
jaccard <- function(a,b) {
  I <- sum((a==1 & b==1)*1)
  U <- sum((a==1 | b==1)*1)
  J <- I/U
  return(J)
}

#### Parameters of bipartite network ####
agents <- 100
artifacts <- 100
density <- .1

#### Compute backbone similarity in bipartites with different degree sequence distributions ####
for (i in 1:100) {  #Number of trials per condition
  for (r_loop in c("left", "right", "normal", "constant", "uniform")) {  #Loop over types of row degree distributions
    if (r_loop=="right") {
      r_alpha <- 1
      r_beta <- 10
    }
    if (r_loop=="left") {
      r_alpha <- 10
      r_beta <- 1
    }
    if (r_loop=="normal") {
      r_alpha <- 10
      r_beta <- 10
    }
    if (r_loop=="constant") {
      r_alpha <- 10000
      r_beta <- 10000
    }
    if (r_loop=="uniform") {
      r_alpha <- 1
      r_beta <- 1
    }
    for (c_loop in c("left", "right", "normal", "constant", "uniform")) {  #Loop over types of column degree distributions
     if (c_loop=="right") {
       c_alpha <- 1
       c_beta <- 10
     }
      if (c_loop=="left") {
        c_alpha <- 10
        c_beta <- 1
      }
      if (c_loop=="normal") {
        c_alpha <- 10
        c_beta <- 10
      }
      if (c_loop=="constant") {
        c_alpha <- 10000
        c_beta <- 10000
      }
      if (c_loop=="uniform") {
        c_alpha <- 1
        c_beta <- 1
      }
      print(paste0("Loop ", i, "   rows ", r_loop, "   columns ", c_loop))
      
      #Create the bipartite
      B <- bipartite.from.distribution(agents,artifacts,density,rowdist=c(r_alpha,r_beta),coldist=c(c_alpha,c_beta))

      #Extract backbones for alpha = 0.05
      ffm.bb <- backbone.extract(fixedfill(B), signed = FALSE, alpha = 0.05)
      frm.bb <- backbone.extract(fixedrow(B), signed = FALSE, alpha = 0.05)
      fcm.bb <- backbone.extract(fixedcol(B), signed = FALSE, alpha = 0.05)
      fdsm.bb <- backbone.extract(fdsm(B, trials = 1000), signed = FALSE, alpha = 0.05)
      sdsm.probs <- sdsm(B)
      sdsm05.bb <- backbone.extract(sdsm.probs, signed = FALSE, alpha = 0.05)
      
      #Extract optimal-alpha SDSM backbone
      for (alpha in seq(from = 0.05, to = .5, by = 0.001)) {
        backbone <- backbone.extract(sdsm.probs, signed = FALSE, alpha = alpha)
        j <- jaccard(backbone,fdsm.bb)
        if (exists("results") == FALSE) {
          results <- data.frame(alpha = alpha, jaccard = j)
        } else results <- rbind(results, c(alpha, j))
      }
      optimal.alpha <- min(results[1][which(results$jaccard==max(results$jaccard)),])
      rm(results)
      sdsmoptim.bb <- backbone.extract(sdsm.probs, signed = FALSE, alpha = optimal.alpha)
      
      #Compute Jaccard
      ffm.j <- jaccard(ffm.bb,fdsm.bb)
      frm.j <- jaccard(frm.bb,fdsm.bb)
      fcm.j <- jaccard(fcm.bb,fdsm.bb)
      sdsm05.j <- jaccard(sdsm05.bb,fdsm.bb)
      sdsmoptim.j <- jaccard(sdsmoptim.bb,fdsm.bb)

      #Count edges
      ffm.e <- sum(ffm.bb)
      frm.e <- sum(frm.bb)
      fcm.e <- sum(fcm.bb)
      sdsm05.e <- sum(sdsm05.bb)
      sdsmoptim.e <- sum(sdsmoptim.bb)
      
      #Post to results object
      if (exists("jaccards") == FALSE) {
        jaccards <- data.frame(loop = i,
                               agents = agents,
                               artifacts = artifacts,
                               density = density,
                               row = r_loop,
                               col = c_loop,
                               ffm.j = ffm.j,
                               ffm.e = ffm.e,
                               frm.j = frm.j,
                               frm.e = frm.e,
                               fcm.j = fcm.j,
                               fcm.e = fcm.e,
                               sdsm05.j = sdsm05.j,
                               sdsm05.e = sdsm05.e,
                               sdsm.alpha = optimal.alpha,
                               sdsmoptim.j = sdsmoptim.j,
                               sdsmoptim.e = sdsmoptim.e)
      } else jaccards <- rbind(jaccards, c(i, agents, artifacts, density, r_loop, c_loop,
                                           ffm.j, ffm.e,
                                           frm.j, frm.e,
                                           fcm.j, fcm.e,
                                           sdsm05.j, sdsm05.e,
                                           optimal.alpha, sdsmoptim.j, sdsmoptim.e))
    }
  }
}
write.csv(jaccards, "study3_jaccards.csv")

#### Compute values to report in manuscript ####
rm(list=ls())
jaccards <- read.csv("study3_jaccards.csv")

m <- mean(jaccards$sdsm.alpha[which(jaccards$row=="constant")])
sd <- sd(jaccards$sdsm.alpha[which(jaccards$row=="constant")])
print(paste0("Mean (sd) alpha to use when rows are constant: ", round(m,3), " (", round(sd,3),")"))

m <- mean(jaccards$sdsm.alpha[which(jaccards$col=="constant")])
sd <- sd(jaccards$sdsm.alpha[which(jaccards$col=="constant")])
print(paste0("Mean (sd) alpha to use when columns are constant: ", round(m,3), " (", round(sd,3),")"))

m <- mean(jaccards$sdsm.alpha[which(jaccards$row!="constant" & jaccards$col!="constant")])
sd <- sd(jaccards$sdsm.alpha[which(jaccards$row!="constant" & jaccards$col!="constant")])
print(paste0("Mean (sd) alpha to use when rows and columns are NOT constant: ", round(m,3), " (", round(sd,3),")"))

jaccards <- jaccards[,c("row", "col", "ffm.j", "frm.j", "fcm.j", "sdsm05.j", "sdsmoptim.j")]
jaccards <- reshape2::melt(jaccards)

j <- mean(jaccards$value[which(jaccards$row=="constant" & (jaccards$col=="constant" | jaccards$row=="left"))])
print(paste0("Mean J when row = cons, col = cons/left: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="frm.j" & jaccards$col=="right")])
print(paste0("Mean J in FRM when col = right: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="frm.j" & jaccards$col=="normal")])
print(paste0("Mean J in FRM when col = normal: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="fcm.j" & (jaccards$row=="right" | jaccards$row=="uniform"))])
print(paste0("Mean J in FCM when row = right or uniform: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="fcm.j" & (jaccards$row=="left" | jaccards$row=="left"))])
print(paste0("Mean J in FCM when row = left or constant: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="fcm.j" & jaccards$col=="uniform")])
print(paste0("Mean J in FCM when col = uniform: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="ffm.j" & (jaccards$col=="right" | jaccards$row=="right" | jaccards$row=="uniform"))])
print(paste0("Mean J in FCM when row = right and/or col = right or uniform: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="ffm.j" & !(jaccards$col=="right" | jaccards$row=="right" | jaccards$row=="uniform"))])
print(paste0("Mean J in FCM when row = right and/or col = right or uniform: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="sdsm05.j" & jaccards$row=="constant" & (jaccards$col=="constant" | jaccards$row=="left"))])
print(paste0("Mean J in SDSM when row = uniform and col = constant or left: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="sdsm05.j" & !(jaccards$row=="constant" & (jaccards$col=="constant" | jaccards$row=="left")))])
print(paste0("Mean J in SDSM when row = not uniform and col = not constant or left: ", round(j,3)))

j <- mean(jaccards$value[which(jaccards$variable=="sdsmoptim.j")])
print(paste0("Mean J in SDSM alpha = optimal: ", round(j,3)))

#### Plotting results for alpha = 0.05 ####
rm(list=ls())
library(ggplot2)
library(cowplot)
library(latex2exp)

# Compute jaccard means within condition
jaccards <- read.csv("study3_jaccards.csv")
jaccards <- aggregate(cbind(ffm.j, frm.j, fcm.j, sdsm05.j, sdsm.alpha, sdsmoptim.j) ~ row + col, data = jaccards, mean)
jaccards$row <- factor(jaccards$row, levels = c("right", "uniform", "constant", "left", "normal"), labels = c("Right", "Unif", "Cons", "Left", "Norm"))
jaccards$col <- factor(jaccards$col, levels = c("right", "uniform", "constant", "left", "normal"), labels = c("Right", "Unif", "Cons", "Left", "Norm"))

# Heatmaps
fixedfill.heatmap <- ggplot(jaccards, aes(y = row, x = col, fill = ffm.j)) + geom_tile() + 
  xlab("Artifact degree distribution") + ylab("Agent degree distribution") + 
  ggtitle("Fixed Fill Model") +
  scale_fill_continuous(limits=c(0, 1), breaks=seq(0,1,by=0.25), guide = guide_legend(override.aes = list(alpha = 0))) + 
  theme(plot.title = element_text(size = 10), axis.title=element_text(size=9), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(angle=90, hjust=.5))

fixedrow.heatmap <- ggplot(jaccards, aes(y = row, x = col, fill = frm.j)) + geom_tile() + 
  xlab("Artifact degree distribution") + ylab("Agent degree distribution") + 
  ggtitle("Fixed Row Model") +
  scale_fill_continuous(limits=c(0, 1), breaks=seq(0,1,by=0.25), guide = guide_legend(override.aes = list(alpha = 0))) + 
  theme(plot.title = element_text(size = 10), axis.title=element_text(size=9), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(angle=90, hjust=.5))

fixedcol.heatmap <- ggplot(jaccards, aes(y = row, x = col, fill = fcm.j)) + geom_tile() + 
  xlab("Artifact degree distribution") + ylab("Agent degree distribution") + 
  ggtitle("Fixed Column Model") +
  scale_fill_continuous(limits=c(0, 1), breaks=seq(0,1,by=0.25), guide = guide_legend(override.aes = list(alpha = 0))) + 
  theme(plot.title = element_text(size = 10), axis.title=element_text(size=9), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(angle=90, hjust=.5))

sdsm05.heatmap <- ggplot(jaccards, aes(y = row, x = col, fill = sdsm05.j)) + geom_tile() + 
  xlab("Artifact degree distribution") + ylab("Agent degree distribution") + labs(fill="Jaccard") + 
  ggtitle("Stochastic Degree Sequence Model") +
  scale_fill_continuous(limits=c(0, 1), breaks=seq(0,1,by=0.25)) + 
  theme(plot.title = element_text(size = 10), axis.title=element_text(size=9), legend.title=element_text(size=9),                                                                       
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(angle=90, hjust=.5))

#Combine heatmaps
heatmaps <- plot_grid(fixedfill.heatmap + theme(legend.position="none"), 
                      fixedrow.heatmap + theme(legend.position="none"),
                      fixedcol.heatmap + theme(legend.position="none"),
                      sdsm05.heatmap + theme(legend.position="none"),
                      nrow=2, labels="AUTO")
legend <- get_legend(sdsm05.heatmap + 
                     theme(legend.key.size = unit(1, 'cm'), 
                           legend.title = element_text(size=15),
                           legend.text = element_text(size=12)))
plot_grid(heatmaps, legend, rel_widths = c(3, .4))
ggsave("study3_plot.pdf", width = 7.5, height = 6, device = "pdf")

#### Supplemental plot of optimal alphas for SDSM backbone ####
optimalalpha.heatmap <- ggplot(jaccards, aes(y = row, x = col, fill = sdsm.alpha)) + geom_tile() + 
  xlab("Artifact degree distribution") + ylab("Agent degree distribution") + labs(fill=TeX("Optimal $\\alpha$")) +
  ggtitle(TeX("SDSM $\\alpha$ to maximize simlarity with FDSM ($\\alpha = 0.05$)")) +
  scale_fill_continuous(limits=c(0.05, .15), breaks=seq(0.05,.15,by=0.025), low="black", high="purple") + 
  theme(plot.title = element_text(size = 9), axis.title=element_text(size=9), legend.title=element_text(size=9),                                                                       
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(angle=90, hjust=.5))

sdsmoptim.heatmap <- ggplot(jaccards, aes(y = row, x = col, fill = sdsmoptim.j)) + geom_tile() + 
  xlab("Artifact degree distribution") + ylab("Agent degree distribution") + labs(fill="Jaccard") + 
  ggtitle(TeX("Similarity of SDSM ($\\alpha = optimal$) and FDSM ($\\alpha = 0.05$)")) +
  scale_fill_continuous(limits=c(0, 1), breaks=seq(0,1,by=0.25)) + 
  theme(plot.title = element_text(size = 10), axis.title=element_text(size=9), legend.title=element_text(size=9),                                                                       
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(angle=90, hjust=.5))

#Combine
plot_grid(sdsmoptim.heatmap, optimalalpha.heatmap, nrow=1, labels="AUTO")
ggsave("study3_optim.pdf", width = 8, height = 3, device = "pdf")

