#### Replication Materials for Rachel Domagalski's Doctoral Thesis ####
#### Written by Rachel Domagalski, Zachary Neal, Bruce Sagan ####
#### Michigan State University ####
#### Chapter 8: Study 2 ####

#### Load packages and set seed ####
rm(list=ls())
library(backbone)
set.seed(5)

#### Example FDSM and SDSM edge weight distributions ####
## Define jaccard coefficient
jaccard <- function(a,b) {
  I <- sum((a==1 & b==1)*1)
  U <- sum((a==1 | b==1)*1)
  J <- I/U
  return(J)
}

## Import GaWC data
gawc <- read.csv(file="da11.csv", header = TRUE, row.names = 1)
gawc <- as.matrix(gawc)
gawc[gawc <= 2] <- 0  #Recode 0s, 1s, and 2s (i.e., Typical and smaller offices) as 0
gawc[gawc >= 3] <- 1  #Recode 3s, 4s, and 5s (i.e., Larger than typical offices) as 1
gawc <- gawc[rowSums(gawc) != 0,]  #Remove rows that are all 0s (cities with no offices)

## Compute edge weight distributions
bicm.prob <- bicm(gawc)
row.names(bicm.prob) <- row.names(gawc)
bicm.weight <- PoissonBinomial::rpbinom(10000, (bicm.prob["MILAN",] * bicm.prob["PARIS",]))
fdsm <- fdsm(gawc, trials = 10000, dyad=c("MILAN","PARIS"))
fdsm.weight <- fdsm$dyad_values
example <- data.frame(null = c(bicm.weight, fdsm.weight), Ensemble = c(rep(c("SDSM"), each = 10000), rep(c("FDSM"), each = 10000)))

## Backbone for SDSM and FDSM
fdsm.bb <- backbone.extract(fdsm, alpha = 0.05)
sdsm <- sdsm(gawc)
sdsm.bb05 <- backbone.extract(sdsm, alpha = 0.05)
sdsm.bb115 <- backbone.extract(sdsm, alpha = 0.115)

## Add supplementary values
R <- nrow(gawc)
C <- ncol(gawc)
N <- sum(gawc)
D <- mean(gawc)
observed <- (gawc%*%t(gawc))["MILAN","PARIS"]
fdsm.p <- as.matrix(fdsm$positive)["MILAN","PARIS"]
sdsm.p <- sdsm$positive["MILAN","PARIS"]
jaccard05 <- jaccard(fdsm.bb,sdsm.bb05)
jaccard115 <- jaccard(fdsm.bb, sdsm.bb115)
fdsm.density <- mean(fdsm.bb)
sdsm05.density <- mean(sdsm.bb05)
sdsm115.density <- mean(sdsm.bb115)
example$values <- c(R, C, N, D, observed, fdsm.p, sdsm.p, jaccard05, jaccard115, fdsm.density, sdsm05.density, sdsm115.density, rep("NA", each = 9988))
write.csv(example,"study2_example.csv")

#### Compare FDSM (alpha = 0.05) backbone to SDSM backbones in GaWC network ####
fdsm.bb <- backbone.extract(fdsm, alpha = 0.05)
sdsm.probs <- sdsm(gawc)

##Extract SDSM at different alpha, comparing each to FDSM
for (alpha in seq(from = 0.01, to = 0.3, by = 0.001)) {
  backbone <- backbone.extract(sdsm.probs, signed = FALSE, alpha = alpha)
  j <- jaccard(backbone,fdsm.bb)
  if (exists("results") == FALSE) {
    results <- data.frame(condition = 0, alpha = alpha,jaccard = j) 
  } else results <- rbind(results, c(0, alpha, j))
}

#### Compare FDSM (alpha = 0.05) backbone to SDSM backbones in simulated GaWC-like networks ####
for (loop in 1:100) {
  print(paste0("Loop ", loop))
  B <- bipartite.from.distribution(196,100,.08331633,rowdist=c(.3,10),coldist=c(.8,10))
  fdsm.bb <- backbone.extract(fdsm(B, trials = 1000), signed = FALSE)
  sdsm.probs <- sdsm(B)

  ##Extract SDSM at different alpha, comparing each to FDSM
  for (alpha in seq(from = 0.01, to = 0.3, by = 0.001)) {
    backbone <- backbone.extract(sdsm.probs, signed = FALSE, alpha = alpha)
    j <- jaccard(backbone,fdsm.bb)
    if (exists("results") == FALSE) {
      results <- data.frame(condition = 1, alpha = alpha,jaccard = j) 
    } else results <- rbind(results, c(1, alpha, j))
  }
}
  
write.csv(results,"study2_power.csv")

#### Plot ####
rm(list=ls())
library(ggplot2)
library(cowplot)
library(latex2exp)

## Distribution example
example <- read.csv("study2_example.csv", row.names = 1, header = TRUE)
observed <- example$values[5]

sink(file="study2_output.txt")
print("Details from Example")
print(paste0("Agents = ", example$values[1]))
print(paste0("Artifacts = ", example$values[2]))
print(paste0("Ones = ", example$values[3]))
print(paste0("Density = ", example$values[4]))
print(paste0("Observed = ", example$values[5]))
print(paste0("FDSM p = ", example$values[6]))
print(paste0("SDSM p = ", example$values[7]))
print(paste0("J(alpha = 0.05) = ", example$values[8]))
print(paste0("J(alpha = 0.115) = ", example$values[9]))
print(paste0("FDSM (alpha = 0.05) Density = ", example$values[10]))
print(paste0("SDSM (alpha = 0.05) Density = ", example$values[11]))
print(paste0("SDSM (alpha = 0.115) Density = ", example$values[12]))
sink()

density.plot <-
ggplot(example, aes(x = null, fill = Ensemble, color = Ensemble)) + 
  geom_density(bw = 1, alpha = .6) +
  geom_vline(xintercept = observed, linetype="dashed", color = "black", size=.5) +
  xlab(TeX("Firms co-located in Milan & Paris (P^*_{ij})")) + ylab("Probability") +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  scale_x_continuous(breaks=c(10,20,30,observed), labels = c("10", "20", "30", TeX("P_{ij} = 26"))) + 
  scale_color_manual(values = c(rgb(1,1,1,1), rgb(1,1,1,1))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.position = c(0.2, 0.8), legend.text=element_text(size=15),
        legend.title=element_text(size=15), legend.key=element_blank())

## Power analysis
results <- read.csv("study2_power.csv")
results <- results[complete.cases(results), ]
observed <- results[which(results$condition==0),][,3:4]

simulated <- results[which(results$condition==1),][,3:4]
lowpct <- function(x) {return(quantile(x,.1))}
highpct <- function(x) {return(quantile(x,.9))}
mean <- aggregate(simulated, by = list(simulated$alpha), mean)
upper <- aggregate(simulated, by = list(simulated$alpha), highpct)[,2:3]
colnames(upper) <- c("alpha", "upper")
lower <- aggregate(simulated, by = list(simulated$alpha), lowpct)[,2:3]
colnames(lower) <- c("alpha", "lower")
ribbon <- merge(upper, lower)

observed.maxjaccard <- max(observed$jaccard)
observed.maxalpha <- min(observed$alpha[which(observed$jaccard==observed.maxjaccard)])
simulated.maxjaccard <- max(mean$jaccard)
simulated.maxalpha <- min(mean$alpha[which(mean$jaccard==simulated.maxjaccard)])

power.plot <-
  ggplot() + 
  geom_line(data = observed, mapping = aes(x = alpha, y = jaccard), 
            color = rgb(124,174,0,maxColorValue = 255), linetype = "solid", size=1.5, alpha = 1) +
  geom_line(data = mean, mapping = aes(x = alpha, y = jaccard), 
            color = rgb(199,124,255,maxColorValue = 255), linetype = "solid", size=1.5, alpha = 1) +
  geom_ribbon(data = ribbon, mapping = aes(x = alpha, ymin = lower, ymax = upper),
              color = NA, fill = rgb(199,124,255,maxColorValue = 255), alpha = 0.5) +
  geom_segment(aes(x = observed.maxalpha, y = 0, xend = observed.maxalpha, yend = observed.maxjaccard), linetype = "dashed") +
  geom_segment(aes(x = 0, y = observed.maxjaccard, xend = observed.maxalpha, yend = observed.maxjaccard), linetype = "dashed") +
  geom_segment(aes(x = simulated.maxalpha, y = 0, xend = simulated.maxalpha, yend = simulated.maxjaccard), linetype = "dashed") +
  geom_segment(aes(x = 0, y = simulated.maxjaccard, xend = simulated.maxalpha, yend = simulated.maxjaccard), linetype = "dashed") +
  expand_limits(y = c(0, 1)) + labs(color="Size, density & degree dist. of" ~ bold("B")) +
  xlab(expression(paste("SDSM Significance level (", alpha, ")"))) + 
  scale_x_continuous(breaks=c(0.05, 0.2, 0.25, 0.3, observed.maxalpha, simulated.maxalpha),
                     labels=c("0.05", "0.20", "0.25", "0.30", round(observed.maxalpha,2), round(simulated.maxalpha,2))) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.00, observed.maxjaccard, simulated.maxjaccard),
                     labels=c("0", "0.25", "0.50", "0.75", "1.00", round(observed.maxjaccard,2), round(simulated.maxjaccard,2))) +
  ylab(TeX("Jaccard simiilarity of $\\mathbf{P'}^{SDSM}_{\\alpha}$ and $\\mathbf{P'}^{FDSM}_{\\alpha = 0.05}$")) +
  annotate("text", x = .045, y = 1, label = "Empirical GaWC network", size = 5, hjust = 0, vjust = .5) +
  annotate("segment", x = 0, xend = .04, y = 1, yend = 1, size = 1.5, colour = rgb(124,174,0,maxColorValue = 255)) +
  annotate("text", x = .045, y = .92, label = "100 simulated networks", size = 5, hjust = 0, vjust = .5) +
  annotate("segment", x = 0, xend = .04, y = .92, yend = .92, size = 1.5, colour = rgb(199,124,255,maxColorValue = 255)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.position = c(0.65, 0.14), legend.text=element_text(size=10),
        legend.title=element_text(size=10), legend.key=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

## Combine
plot_grid(density.plot, power.plot, ncol=2, labels="AUTO")
ggsave("study2_plot.pdf", width = 12, height = 6, device = "pdf")
