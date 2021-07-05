#### Replication Materials for Rachel Domagalski's Doctoral Thesis ####
#### Written by Rachel Domagalski, Zachary Neal, Bruce Sagan ####
#### Michigan State University ####
#### Chapter 8: Study 1 ####

rm(list=ls())
library(backbone)
set.seed(5)

##############################
### BEGIN HELPER FUNCTIONS ###
##############################

#### Compute Pr(B_ij=1) in all B with fixed row/column degree sequences ####
true_probs <- function(B, progress = FALSE){  
  upperbound.row <- 1  #Find row-wise upper bound on cardinality
  for (i in 1:nrow(B)) {
    thisrow <- choose(ncol(B), sum(B[i,]))  #Ways to place sum(B[i,]) 1s into ncol(B) columns
    upperbound.row <- upperbound.row * thisrow
  }
  upperbound.col <- 1  #Find column-wise upper bound on cardinality
  for (i in 1:ncol(B)) {
    thiscol <- choose(nrow(B), sum(B[,i]))  #Ways to place sum(B[,i]) 1s into nrow(B) rows
    upperbound.col <- upperbound.col * thiscol
  }
  upperbound <- min(upperbound.row, upperbound.col)  #A generous upper bound on cardinality
  trials <- (upperbound * 5) - 1  #Generate 5 times this many random B (minus 1, because we can also use B)
  Bstar <- t(as.vector(B))  #Row vector containing B
  
  if (progress == TRUE) {pb <- txtProgressBar(min = 0, max = trials, style = 3)}
  for (i in 1:trials) {  #Add a new row vector for each random B generated using Curveball algorithm
    Bstar <- rbind(Bstar,t(as.vector(backbone::curveball(B))))
    Bstar <- unique(Bstar)  #Only keep a new row if it's unique
    if (progress == TRUE) {setTxtProgressBar(pb, i)}
  }
  
  probs <- colSums(Bstar)/nrow(Bstar)  #Compute probability that B_ij = 1 in this space
  probs <- matrix(probs, nrow = nrow(B), ncol = ncol(B))  #Reassemble as matrix
  return(list("probability" = probs, "cardinality" = nrow(Bstar)))
} #end true_probs

#### Helper LDM function from https://statisticalhorizons.com/better-predicted-probabilities ####
predict_ldm <- function(fit) {
  c <- fit$coefficients[1]
  k <- nrow( fit$model ) / deviance( fit )
  m <- mean( fit$model[,1] )
  a <- log(m/(1-m) ) + k*( c-.5 ) + .5*(1/m - 1/(1-m) )
  lin_preds <- predict(fit)
  my_preds_logit <- k*(lin_preds-c) + a
  my_preds <- 1/(1 + exp(-my_preds_logit))
  names(my_preds) <- "ldm_preds"
  return(my_preds)
}

#### Polytope model ####
polytope <- function(G){
  
  #### Define Variable to solve for
  matrix <- CVXR::Variable(dim(G)[1], dim(G)[2])
  
  #### Define row & column sums
  mat1 <- matrix(1, dim(G)[2], 1)
  mat2 <- matrix(1, dim(G)[1], 1)
  
  #### Define Constraints
  constraint1 <- matrix >= 0
  constraint2 <- matrix <= 1
  if (methods::is(G, "sparseMatrix")) {constraint3 <- (matrix%*%mat1) == Matrix::rowSums(G)
  } else {constraint3 <- (matrix%*%mat1) == rowSums(G)}
  if (methods::is(G, "sparseMatrix")) {constraint4 <- t(matrix)%*%mat2 == Matrix::colSums(G)
  } else {constraint4 <- t(matrix)%*%mat2 == colSums(G)}
  constraints <- list(constraint1, constraint2, constraint3, constraint4)
  
  #### Define Objective, the function to solve
  objective <- CVXR::Maximize(sum(CVXR::entr(matrix)+CVXR::entr(1-matrix)))
  
  #### Define Problem: objective with the constrants
  problem <- CVXR::Problem(objective, constraints)
  
  #### Solve the problem ####
  result <- suppressWarnings(CVXR::psolve(problem))
  
  ### Warning/Stop if not optimal ###
  if (result$status == "optimal_inaccurate") {warning("polytope result not optimal")}
  if (result$status != "optimal" & result$status != "optimal_inaccurate") {stop("unable to compute SDSM-Polytope")}
  
  #### Results
  new_matrix <- result$getValue(matrix)
  
  ### Restrict values between 0 and 1 ###
  gr <- which(new_matrix>1)
  new_matrix[gr] <- 1
  le <- which(new_matrix<0)
  new_matrix[le] <- 0
  
  #### Return Matrix of Probabilities ####
  return(new_matrix)
}

#### Estimate Pr(B_ij=1) in all B with fixed row/column marginals ####
library(Rcpp)
sourceCpp("scobit.cpp") #Import scobit model estimation, implemented in C++ by David Schoch
sourceCpp("scobit_interaction.cpp") #Import scobit model estimation, implemented in C++ by David Schoch

est_probs <- function(B, model){  
  
  ### Compute row and column sums if necessary ###
  if (model=="logit" | model=="logit.i" | model=="scobit" | model=="scobit.i" | model=="lpm" | model=="ldm" | model=="rcn") {
    ## Vectorize the bipartite data ##
    A <- data.frame(as.vector(B))
    names(A)[names(A)=="as.vector.B."] <- "value"
    
    ## Assign row and column IDs in the vectorized data ##
    A$row <- rep(1:nrow(B), times=ncol(B))
    A$col <- rep(1:ncol(B), each=nrow(B))
    
    ## Compute and attach rowsums, columnsums ##
    A$rowmarg <- stats::ave(A$value,A$row,FUN=sum)
    A$colmarg <- stats::ave(A$value,A$col,FUN=sum)
    A$rowcol <- A$rowmarg * A$colmarg
  }
  
  ### Logit model without interaction ###
  if (model=="logit") {
    model.estimates <- speedglm::speedglm(formula= value ~  rowmarg + colmarg, family = stats::binomial(link="logit"), data=A)
    probs <- stats::predict(model.estimates,newdata=A,type = "response")
    probs <- matrix(probs, nrow=nrow(B), ncol=ncol(B))
  }
  
  ### Logit model with interaction ###
  if (model=="logit.i") {
    model.estimates <- speedglm::speedglm(formula= value ~  rowmarg + colmarg + rowcol, family = stats::binomial(link="logit"), data=A)
    probs <- stats::predict(model.estimates,newdata=A,type = "response")
    probs <- matrix(probs, nrow=nrow(B), ncol=ncol(B))
  }
  
  ### Scobit model without interaction ###
  if (model == "scobit") {
    params <- list(b0=0.1,b1=0.00005,b2=0.00005,a=0.01)
    model.estimates <- stats::optim(params,scobit_loglike_cpp,gr=scobit_loglike_gr_cpp,method="BFGS",x1=A$rowmarg,x2=A$colmarg,y=A$value)
    pars <- c(model.estimates$par[1],model.estimates$par[2],model.estimates$par[3])
    probs <- 1-1/(1+exp(pars[1]+pars[2]*A$rowmarg+pars[3]*A$colmarg))^model.estimates$par[4]
    probs <- matrix(probs, nrow=nrow(B), ncol=ncol(B))
  }
  
  ### Scobit model with interaction ###
  if (model == "scobit.i") {
    params <- list(b0=0.1,b1=0.00005,b2=0.00005, b3 = 0.00005, a=0.01)
    model.estimates <- stats::optim(params,scobit_loglike_cpp_i,gr=scobit_loglike_gr_cpp_i,method="BFGS",x1=A$rowmarg,x2=A$colmarg,y=A$value)
    pars <- c(model.estimates$par[1],model.estimates$par[2],model.estimates$par[3], model.estimates$par[4])
    probs <- 1-1/(1+exp(pars[1]+pars[2]*A$rowmarg+pars[3]*A$colmarg+pars[4]*A$rowcol))^model.estimates$par[5]
    probs <- matrix(probs, nrow=nrow(B), ncol=ncol(B))
  }
  
  ### Linear probability model ###
  if (model=="lpm") {
    model.estimates <- stats::lm(formula= value ~ rowmarg + colmarg, data=A)
    probs <- stats::predict(model.estimates)
    probs[probs<0] <- 0 #Truncate out-of-bounds estimates
    probs[probs>1] <- 1
    probs <- matrix(probs, nrow=nrow(B), ncol=ncol(B))
  }
  
  ### Linear discriminant model ###
  if (model=="ldm") {
    model.estimates <- stats::lm(formula= value ~ rowmarg + colmarg, data=A)
    probs <- predict_ldm(model.estimates)
    probs[probs<0] <- 0 #Truncate out-of-bounds estimates
    probs[probs>1] <- 1
    probs <- matrix(probs, nrow=nrow(B), ncol=ncol(B))
  }
  
  ### Chi-Square model ###
  if (model=="rcn") {
    probs <- (A$rowmarg * A$colmarg)/sum(A$value)
    probs[probs<0] <- 0 #Truncate out-of-bounds estimates
    probs[probs>1] <- 1
    probs <- matrix(probs, nrow=nrow(B), ncol=ncol(B))
  }
  
  ### Polytopes model ###
  if (model=="polytope") {
    probs <- polytope(B)
  }
  
  ### Bipartite Configuation model ###
  if (model=="bicm") {
    probs <- bicm(B)
  }
  
  return(probs)
}
############################
### END HELPER FUNCTIONS ###
############################

#### Identify true Pr(B_ij = 1) for all B_FDSM ####
## All possible degree sequences
samples <- 100000
print("Identifying all possible degree sequences of matrices sizes 3x3 to 5x5")
for (loop in 1:samples) {
  r <- sample(c(3,4,4,4,4,5,5,5,5,5,5,5,5,5,5), 1)  #Sample a number of rows, with greater chance of picking larger number
  if (r == 5) {c <- 5} else {c <- sample(r:5, 1)}  #Sample a number of columns, equal or greater than rows
  matrix <- bipartite.from.probability(r,c)
  R <- sort(rowSums(matrix))  #Compute row sums, sort them
  C <- sort(colSums(matrix))  #Compute column sums, sort them
  if (exists("sequences") == FALSE) {  #Put each set of sorted sequences in the "sequences" object
    sequences <- data.frame(rowsums = toString(R),
                            colsums = toString(C))
  } else sequences <- rbind(sequences, c(rowsums = toString(R),
                                         colsums = toString(C)))
  sequences <- unique(sequences)  #Keep only unique sequences  
  if (loop == 1) {
    print("Identifying all possible degree sequences")
    pb <- txtProgressBar(min = 0, max = samples, style = 3)
  }
  setTxtProgressBar(pb, loop)
}
sequences <- sequences[order(nchar(sequences$rowsums), nchar(sequences$colsums)),]

## Probabilities
for (loop in 1:nrow(sequences)) {
  
  ### Generate a matrix with the prescribed degree sequence
  R <- as.numeric(unlist(strsplit(sequences$rowsums[loop], split=", ")))  #Read in a row degree sequence
  C <- as.numeric(unlist(strsplit(sequences$colsums[loop], split=", ")))  #Read in a column degree sequence
  matrix <- bipartite.from.sequence(R, C)  #Create a bipartite with the requested row and column degree sequences
  
  ###Obtain true probabilities
  true.prob <- true_probs(matrix, progress = TRUE)  #Compute Pr(B_ij = 1) for all B in B(R,C)
  
  ###Post to results object
  if (exists("probabilities") == FALSE) {
    probabilities <- data.frame(rows = length(R),
                                rowsums = toString(R),
                                cols = length(C),
                                colsums = toString(C),
                                cardinality = true.prob$cardinality,
                                probs = toString(as.vector(true.prob$probability)))
  } else probabilities <- rbind(probabilities, c(rows = length(R),
                                                 rowsums = toString(R),
                                                 cols = length(C),
                                                 colsums = toString(C),
                                                 cardinality = true.prob$cardinality,
                                                 probs = toString(as.vector(true.prob$probability))))
  print(loop)
}
write.csv(probabilities, "study1_trueprobs.csv")

#### Accuracy: Compare true and estimated probabilities ####
true_probs <- read.csv("study1_trueprobs.csv", header = TRUE, row.names = 1)

for (loop in 1:nrow(true_probs)) {
  ### Read in a degree sequence and associated true probabilities
  R <- as.numeric(unlist(strsplit(true_probs$rowsums[loop], split=", ")))  #Read in a row degree sequence
  C <- as.numeric(unlist(strsplit(true_probs$colsums[loop], split=", ")))  #Read in a column degree sequence
  card <- true_probs$cardinality[loop]  #Read in the cardinality
  true <- as.numeric(unlist(strsplit(true_probs$probs[loop], split=", ")))  #Read in the true probabilities
  
  ### Generate a matrix with the prescribed degree sequence
  matrix <- bipartite.from.sequence(R, C)
  
  ### Compute estimated probabilities
  logit.prob <- est_probs(matrix, model = "logit")
  scobit.prob <- est_probs(matrix, model = "scobit")
  logit.i.prob <- est_probs(matrix, model = "logit.i")
  scobit.i.prob <- est_probs(matrix, model = "scobit.i")
  lpm.prob <- est_probs(matrix, model = "lpm")
  ldm.prob <- est_probs(matrix, model = "ldm")
  rcn.prob <- est_probs(matrix, model = "rcn")
  poly.prob <- est_probs(matrix, model = "polytope")
  bicm.prob <- est_probs(matrix, model = "bicm")

  ### Post to results object
  if (exists("results") == FALSE) {
  results <- data.frame(rows = length(R),
                        rowsums = toString(R),
                        cols = length(C),
                        colsums = toString(C),
                        cardinality = card,
                        scobit.i = mean(abs(true - scobit.i.prob)),
                        logit.i = mean(abs(true - logit.i.prob)),
                        rcn = mean(abs(true - rcn.prob)),
                        lpm = mean(abs(true - lpm.prob)),
                        scobit = mean(abs(true - scobit.prob)),
                        ldm = mean(abs(true - ldm.prob)),
                        logit = mean(abs(true - logit.prob)),
                        poly = mean(abs(true - poly.prob)),
                        bicm = mean(abs(true - bicm.prob)))
  } else results <- rbind(results, c(rows = length(R),
                                     rowsums = toString(R),
                                     cols = length(C),
                                     colsums = toString(C),
                                     cardinality = card,
                                     scobit.i = mean(abs(true - scobit.i.prob)),
                                     logit.i = mean(abs(true - logit.i.prob)),
                                     rcn = mean(abs(true - rcn.prob)),
                                     lpm = mean(abs(true - lpm.prob)),
                                     scobit = mean(abs(true - scobit.prob)),
                                     ldm = mean(abs(true - ldm.prob)),
                                     logit = mean(abs(true - logit.prob)),
                                     poly = mean(abs(true - poly.prob)),
                                     bicm = mean(abs(true - bicm.prob))))
  print(loop)
}
write.csv(results, "study1_accuracy.csv")

#### Speed: Running time of four most accurate estimation methods ####
for (elements in c(2, 2, 2.25, 2.5, 2.75,  #First one loads all packages, so running time is inflated
                      3, 3.25, 3.5, 3.75,
                      4, 4.25, 4.5, 4.75,
                      5, 5.25, 5.5, 5.75,
                      6, 6.25, 6.5, 6.75, 7)){  

  multiplier <- sample(seq(.5,1.5,.1),1)  #Ratio of rows to columns
  r <- round(sqrt(10^elements) * multiplier)  #Number of rows
  c <- round(10^elements / r)  #Number of columns
  p <- runif(1,.25,.75)  #Density
  B <- matrix(rbinom(r*c,1,p),r,c)  #Generate matrix
  print(paste0("Rows = ", r ,"  Columns = ", c , "  Elements = ", r*c))
  
  ### Compute probabilities
  start <- Sys.time()
  logit.prob <- est_probs(B, model = "logit")
  end <- Sys.time()
  logit.time <- as.numeric(difftime(end, start, units = "secs"))
  
  start <- Sys.time()
  ldm.prob <- est_probs(B, model = "ldm")
  end <- Sys.time()
  ldm.time <- as.numeric(difftime(end, start, units = "secs"))
  
  start <- Sys.time()
  poly.prob <- est_probs(B, model = "polytope")
  end <- Sys.time()
  poly.time <- as.numeric(difftime(end, start, units = "secs"))
  
  start <- Sys.time()
  bicm.prob <- est_probs(B, model = "bicm")
  end <- Sys.time()
  bicm.time <- as.numeric(difftime(end, start, units = "secs"))
  
  ### Post to results object
  if (exists("running.time") == FALSE) {
    running.time <- data.frame(r = r,
                               c = c,
                               p = p,
                               Logit = logit.time,
                               LDM = ldm.time,
                               Poly = poly.time,
                               BiCM = bicm.time)
  } else running.time <- rbind(running.time, c(r, c, p, logit.time, ldm.time, poly.time, bicm.time))
}
write.csv(running.time, "study1_time.csv")

#### Generate Plots ####
rm(list=ls())
library(ggplot2)
library(scales)
library(cowplot)
library(latex2exp)
set.seed(5)

## Accuracy
data <- read.csv("study1_accuracy.csv")
deviations <- data[,c("logit.i", "scobit.i", "rcn", "lpm", "scobit", "ldm", "logit", "bicm", "poly")]
means <- data.frame(x = c(1:9), mean = sapply(deviations,mean))
deviations$id <- as.factor(rownames(deviations))
deviations <- reshape2::melt(deviations, id.vars = c("id"))
deviations$variable <- as.character(deviations$variable)
deviations$variable[which(deviations$variable=="logit.i")] <- 1
deviations$variable[which(deviations$variable=="scobit.i")] <- 2
deviations$variable[which(deviations$variable=="rcn")] <- 3
deviations$variable[which(deviations$variable=="lpm")] <- 4
deviations$variable[which(deviations$variable=="scobit")] <- 5
deviations$variable[which(deviations$variable=="ldm")] <- 6
deviations$variable[which(deviations$variable=="logit")] <- 7
deviations$variable[which(deviations$variable=="poly")] <- 8
deviations$variable[which(deviations$variable=="bicm")] <- 9
deviations$variable <- as.numeric(deviations$variable)

accuracy <- 
  ggplot() +
  geom_line(deviations, mapping = aes(x = variable, y = value, group = id), colour = rgb(0,0,0,.06)) +
  geom_line(means, mapping = aes(x = x, y = mean), colour = rgb(1,0,0,1), size = 1) +
  scale_x_continuous(breaks = c(1:9), labels=c("Logit-I", "Scobit-I", "RCF", "LPM", "Scobit", "LDM", "Logit", "Poly", "BiCM")) +
  scale_y_continuous(trans = "sqrt", breaks = c(.01, .05, .1, .2, .3, .4)) +
  annotate("text", x = 4.5, y = .4, label = TeX("One \\textit{B}^{FDSM}"), size = 5, hjust = 0, vjust = .5) +
  annotate("segment", x = 3.5, xend = 4.2, y = .4, yend = .4, size = 1, colour = rgb(0,0,0,.3)) +
  annotate("text", x = 4.5, y = .35, label = TeX("Mean for all 384 \\textit{B}^{FDSM}"), size = 5, hjust = 0, vjust = .5) +
  annotate("segment", x = 3.5, xend = 4.2, y = .35, yend = .35, size = 1, colour = rgb(1,0,0,1)) +
  xlab("Probability estimation method") + ylab("Accuracy: Mean absolute difference") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12), axis.title=element_text(size=15))

## Speed
data <- read.csv("study1_time.csv")[-1,-1]
data$elements <- data$r * data$c
data <- data[,c("elements", "Logit", "LDM", "Poly", "BiCM")]
time <- reshape2::melt(data, id.vars = c("elements"))
colnames(time) <- c("elements", "Method", "time")

speed <- 
  ggplot(time, aes(x=elements, y=time, color=Method)) +
  geom_point(shape = 1) + 
  geom_smooth(method = loess, se = FALSE, span = 500) +
  scale_y_continuous(trans = "log10", breaks = c(.001, .01, .1, 1, 10, 100, 1000),
                     labels = c(expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3))) + 
  scale_x_continuous(trans = "log10", breaks = c(100, 1000, 10000, 100000, 1000000, 10000000, 100000000),
                     labels = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7), expression(10^8))) + 
  xlab("Number of probabilities to estimate") + ylab(expression("Seconds to estimate probabilities")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.position = c(0.2, 0.8), legend.text=element_text(size=15), legend.title=element_blank(),legend.key=element_blank())

## Combine
plot_grid(accuracy, speed, ncol=2, labels="AUTO")
ggsave("study1_plot.pdf", width = 12, height = 5, device = "pdf")

