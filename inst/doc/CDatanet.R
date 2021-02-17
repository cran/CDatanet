## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)

## ----dataexo, echo = TRUE, eval = TRUE----------------------------------------
set.seed(2020) 
library(CDatanet)
# Groups' size
M      <- 5 # Number of sub-groups
nvec   <- round(runif(M, 100, 1000))
print(nvec)
n      <- sum(nvec)
print(n)

# Parameters
lambda <- 0.4              # peer effects
beta   <- c(2, -1.9, 0.8)  # own effects
gamma  <- c(1.5, -1.2)     # contextual effects
sigma  <- 1.5              # standard deviation of epsilon
theta  <- c(lambda, beta, gamma, sigma)
  
# X
data   <- data.frame(x1 = rnorm(n, 1, 1), x2 =  rexp(n, 0.4))

## ----netexo, echo = TRUE, eval = TRUE-----------------------------------------
# Network
Glist  <- list()
for (m in 1:M) {
  nm           <- nvec[m]
  Gm           <- matrix(0, nm, nm)
  max_d        <- 30
  for (i in 1:nm) {
    tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
    Gm[i, tmp] <- 1
  }
  rs           <- rowSums(Gm); rs[rs == 0] <- 1
  Gm           <- Gm/rs
  Glist[[m]]   <- Gm
}

## ----dataexo2, echo = TRUE, eval = TRUE, fig.height = 3, fig.align = "center"----
ytmp    <- simCDnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist, theta = theta,
                    data = data)
names(ytmp)
y       <- ytmp$y
# Add y to data
data$y  <- y
# Summarize y
summary(y)
table(y)
# Plot data histogram
library(ggplot2)
print(ggplot(data = data.frame(y = y), aes(x = y)) +
        geom_bar(color = "black", fill = "cyan") +
        theme_bw() + xlab("") + ylab("Frequency"))

## ----estexoload, echo = FALSE, eval = TRUE------------------------------------
#  I will not run the estimates because it takes time. I will just load the saved ouptuts. However, you will get the same results if you run the code on your computer with the seed set above.
load("out.est.Rdata")

## ----estexo1, echo = TRUE, eval = FALSE---------------------------------------
#  #  Count data
#  CD   <- CDnetNPL(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
#                    optimizer = "nlm",
#                   data = data, npl.ctr = list(print = FALSE, maxit = 5e3))
#  summary(CD)

## ----estexo1p, echo = FALSE, eval = TRUE--------------------------------------
summary(CD)

## ----estexo2, echo = TRUE, eval = FALSE---------------------------------------
#  # SART
#  SART   <- SARTML(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
#                  optimizer = "nlm", data = data, print = FALSE)
#  summary(SART)

## ----estexo2p, echo = FALSE, eval = TRUE--------------------------------------
summary(SART)

## ----estexo3, echo = TRUE, eval = FALSE---------------------------------------
#  #SAR
#  SAR  <- SARML(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
#                optimizer = "nlm", data = data, print = FALSE)
#  summary(SAR)

## ----estexo3p, echo = FALSE, eval = TRUE--------------------------------------
summary(SAR)

## ----mcfun, echo = TRUE, eval = TRUE------------------------------------------
fMC <- function(s) {
  # Groups' size
  M      <- 5
  nvec   <- rep(500, 5); 
  n      <- sum(nvec)
  # Parameters
  lambda <- 0.4 
  beta   <- c(2, -1.9, 0.8)
  gamma  <- c(1.5, -1.2)
  sigma  <- 1.5   
  theta  <- c(lambda, beta, gamma, sigma)
  
  # X
  data   <- data.frame(x1 = rnorm(n, 1, 1), x2 =  rexp(n, 0.4))
  
  # Network
  Glist  <- list()
  for (m in 1:M) {
    nm           <- nvec[m]
    Gm           <- matrix(0, nm, nm)
    max_d        <- 30
    for (i in 1:nm) {
      tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
      Gm[i, tmp] <- 1
    }
    rs           <- rowSums(Gm); rs[rs == 0] <- 1
    Gm           <- Gm/rs
    Glist[[m]]   <- Gm
  }
  
  # y
  ytmp   <- simCDnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist, theta = theta,
                     data = data)
  y      <- ytmp$y
  data$y <- y
  
  #  Models
  CD     <- CDnetNPL(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
                     optimizer = "nlm", 
                     data = data, npl.ctr = list(print = FALSE, maxit = 5e3))
  SART   <- SARTML(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
                   optimizer = "nlm", data = data, print = FALSE)
  SAR    <- SARML(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
                  optimizer = "nlm", data = data, print = FALSE)
  c(CD$estimate, SART$estimate, SAR$estimate)
}

## ----mc, echo = TRUE, eval = FALSE, message = FALSE---------------------------
#  library(doParallel)
#  n.cores  <- 32
#  replic   <- 1000 #Number of replications
#  out.mc   <- mclapply(1:replic, fMC, mc.cores = n.cores)
#  out.mc   <- apply(t(do.call(cbind, out.mc)), 2, mean)
#  out.mc   <- cbind(theta, out.mc[1:7], out.mc[8:14], out.mc[15:21])
#  colnames(out.mc) <- c("TrueValue", "CoutData", "SART", "SAR")
#  print(out.mc)

## ----mcprint, echo = FALSE, eval = TRUE---------------------------------------
load("out.mc.Rdata")
print(out.mc)

## ----dataendo, echo = TRUE, eval = TRUE---------------------------------------
rm(list = ls())
set.seed(2020) 
# Groups' size
M      <- 3 # Number of sub-groups
nvec   <- round(runif(M, 100, 500))
print(nvec)
n      <- sum(nvec)
print(n)

# Parameters
lambda <- 0.4               # peer effects
beta   <- c(2, -1.9, 0.8)   # own effects
gamma  <- c(1.5, -1.2)      # contextual effects
sigma  <- 1.5               # standard deviation of epsilon
theta  <- c(lambda, beta, gamma, sigma)

# X
data   <- data.frame(x1 = rnorm(n, 1, 1), x2 =  rexp(n, 0.4))

## ----netendo, echo = TRUE, eval = TRUE----------------------------------------
# Parameter for network model
betanet  <- c(-2.8, -1.5)    # beta
Glist    <- list()           # adjacency matrix row normalized
Network  <- list()           # adjacency matrix row non-normalized
dX       <- matrix(0, 0, 2)  # observed dyad-specific variables
mu       <- list()           # unobserved individual-level attribute
uu       <- runif(M, -1, 1)  # mean of uu in each sub-network
sigma2u  <- runif(M, 0.5, 4) # variance of uu in each sub-network

# Network
for (m in 1:M) {
  nm           <- nvec[m]
  mum          <- rnorm(nm, uu[m], sqrt(sigma2u[m]))
  Z1           <- matrix(0, nm, nm)  
  Z2           <- matrix(0, nm, nm)
  
  for (i in 1:nm) {
    for (j in 1:nm) {
      Z1[i, j] <- abs(data$x1[i] - data$x1[j])
      Z2[i, j] <- abs(data$x2[i] - data$x2[j])
    }
  }
  
  Gm           <- 1*((Z1*betanet[1] + Z2*betanet[2] +
                        kronecker(mum, t(mum), "+") + rlogis(nm^2)) > 0)
  diag(Gm)     <- 0
  
  diag(Z1)     <- NA
  diag(Z2)     <- NA
  Z1           <- Z1[!is.na(Z1)]
  Z2           <- Z2[!is.na(Z2)]
  
  dX           <- rbind(dX, cbind(Z1, Z2))
  
  Network[[m]] <- Gm
  rs           <- rowSums(Gm); rs[rs == 0] <- 1
  Gm           <- Gm/rs
  Glist[[m]]   <- Gm
  mu[[m]]      <- mum
}
mu            <- unlist(mu)

## ----dataendo2, echo = TRUE, eval = TRUE, fig.height = 3, fig.align = "center"----
tmu      <- (mu - rep(uu, nvec))/sqrt(rep(sigma2u, nvec))
data$tmu <- tmu
rho      <- 0.24
rhobar   <- 0.18
thetanet <- c(lambda, beta, sigma*rho, gamma, sigma*rhobar,
              sigma*(1 - rho - rhobar))
ytmp     <- simCDnet(formula = ~ x1 + x2 + tmu | x1 + x2 + tmu,
                     Glist = Glist, theta = thetanet, data = data)
y        <- ytmp$y
# Add y to data
data$y   <- y
# Summarize y
summary(y)
# Plot data histogram
print(ggplot(data = data.frame(y = y), aes(x = y)) +
        geom_bar(color = "black", fill = "cyan") +
        theme_bw() + xlab("") + ylab("Frequency"))

## ----estendo1, echo = TRUE, eval = TRUE---------------------------------------
# Count data model
CDexo <- CDnetNPL(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
                  optimizer = "optim", data = data, npl.ctr = list(print = FALSE))
summary(CDexo)

## ----estendo2, echo = TRUE, eval = FALSE--------------------------------------
#  # Dyadic linking model
#  net <- netformation(network =  Network, formula = ~ dX, fixed.effects = TRUE,
#                      mcmc.ctr = list(burnin = 1000, iteration = 5000))

## ----mcprintprime, echo = FALSE, eval = TRUE----------------------------------
load("out.net.Rdata")
# I copied and pasted the output to save time. But if you run the code on your laptop with the same seed, you will get the same results
cat("0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|


The program successfully executed 

********SUMMARY******** 
n.obs          :  314890 
n.links        :  26325 
K              :  2 
Fixed effects  :  Yes 
Burnin         :  1000 
Iteration      :  6000 

Elapsed time   :  0  HH  7  mm  56  ss 
 
Average acceptance rate 
                   beta:  0.274125 
                     mu:  0.2697619 ")

## ----estendo2a, echo = TRUE, eval = FALSE-------------------------------------
#  # plot simulations
#  par(mfrow = c(4,2), mar = c(2, 2, 2, 1.9))
#  plot(net$posterior$beta[,1], type = "l", ylim = c(-2.9, -2.6), col = "blue",
#       main = bquote(beta[1]), ylab = "")
#  abline(h = betanet[1], col = "red")
#  
#  plot(net$posterior$beta[,2], type = "l", ylim = c(-1.6, -1.4), col = "blue",
#       main = bquote(beta[2]), ylab = "")
#  abline(h = betanet[2], col = "red")
#  
#  plot(net$posterior$mu[,10], type = "l", col = "blue",
#       main = bquote(mu[10]), ylab = "")
#  abline(h = mu[10], col = "red")
#  
#  plot(net$posterior$mu[,542], type = "l", col = "blue",
#       main = bquote(mu[542]), ylab = "")
#  abline(h = mu[542], col = "red")
#  
#  plot(net$posterior$mu[,849], type = "l", col = "blue",
#       main = bquote(mu[849]), ylab = "")
#  abline(h = mu[849], col = "red")
#  
#  plot(net$posterior$mu[,752], type = "l", col = "blue",
#       main = bquote(mu[752]), ylab = "")
#  abline(h = mu[752], col = "red")
#  
#  plot(net$posterior$uu[,1], type = "l", col = "blue",
#       main = bquote(u[mu][1]), ylab = "")
#  abline(h = uu[1], col = "red")
#  
#  plot(net$posterior$sigmamu2[,3], type = "l", col = "blue",
#       main =  bquote(sigma[mu][3]^2), ylab = "")
#  abline(h = sigma2u[3], col = "red")

## ----estendo2aa, echo = FALSE, eval = TRUE------------------------------------
# plot simulations
par(mfrow = c(4,2), mar = c(2, 2, 2, 1.9))
plot(snet[[1]], type = "l", ylim = c(-2.9, -2.6), col = "blue",
     main = bquote(beta[1]), ylab = "")
abline(h = betanet[1], col = "red")

plot(snet[[2]], type = "l", ylim = c(-1.6, -1.4), col = "blue",
     main = bquote(beta[2]), ylab = "")
abline(h = betanet[2], col = "red")

plot(snet[[3]], type = "l", col = "blue",
     main = bquote(mu[10]), ylab = "")
abline(h = mu[10], col = "red")

plot(snet[[4]], type = "l", col = "blue",
     main = bquote(mu[542]), ylab = "")
abline(h = mu[542], col = "red")

plot(snet[[5]], type = "l", col = "blue",
     main = bquote(mu[849]), ylab = "")
abline(h = mu[849], col = "red")

plot(snet[[6]], type = "l", col = "blue",
     main = bquote(mu[752]), ylab = "")
abline(h = mu[752], col = "red")

plot(snet[[7]], type = "l", col = "blue",
     main = bquote(u[mu][1]), ylab = "")
abline(h = uu[1], col = "red")

plot(snet[[8]], type = "l", col = "blue",
     main =  bquote(sigma[mu][3]^2), ylab = "")
abline(h = sigma2u[3], col = "red")

## ----estendo3, echo = TRUE, eval = FALSE--------------------------------------
#  t           <- which.max(net$posterior$log.density)
#  print(t)

## ----estendo3a, echo = FALSE, eval = TRUE-------------------------------------
print(t)

## ----estendo31, echo = TRUE, eval = FALSE-------------------------------------
#  muest       <- net$posterior$mu[t,]
#  uuest       <- net$posterior$uu[t,]
#  sigma2uest  <- net$posterior$sigmamu2[t,]
#  tmuest      <- (muest - rep(uuest, nvec))/sqrt(rep(sigma2uest, nvec))
#  data$tmuest <- tmuest
#  CDendo      <- CDnetNPL(formula = y ~ x1 + x2 + tmuest, contextual = TRUE,
#                          Glist = Glist, optimizer = "optim", data = data,
#                          npl.ctr = list(print = FALSE))
#  summary(CDendo)

## ----estendo31a, echo = FALSE, eval = TRUE------------------------------------
data$tmuest <- tmuest
summary(CDendo)

## ----estendo4, echo = TRUE, eval = FALSE--------------------------------------
#  fendo    <- function(s) {
#    t         <- sample(3001:8000, 1)
#    datat     <- data
#    mus       <- net$posterior$mu[t,]
#    uus       <- rep(net$posterior$uu[t,], nvec)
#    sus       <- rep(net$posterior$sigmamu2[t,], nvec)
#    tmu       <- (mus - uus)/sqrt(sus)
#    datat$tmu <- tmu
#  
#    CDnet    <- CDnetNPL(formula = y ~ x1 + x2 + tmu, contextual = TRUE,
#                         Glist = Glist, optimizer = "optim", data = datat,
#                         npl.ctr = list(print = FALSE))
#  
#    summary(CDnet, Glist = Glist, data = datat)
#  }
#  
#  n.cores    <- 32
#  replic     <- 1000 #Number of replications
#  out.endo   <- mclapply(1:replic, fendo, mc.cores = n.cores)
#  # The output of out.endo is a list of objects from "summary.CDnetNPL" class
#  # Let's set the class of out.endo as "summary.CDnetNPLs"
#  class(out.endo) <- "summary.CDnetNPLs"
#  # I can now summarize the results using print
#  # Object from "summary.CDnetNPL" has a print method
#  print(out.endo)

## ----endoprint1, echo = FALSE, eval = TRUE------------------------------------
load("out.endo.Rdata")
CDatanet:::print.summary.CDnetNPLs(out.endo)

## ----mcfunendo, echo = TRUE, eval = TRUE--------------------------------------
fMCendo <- function(s) {
  # parameters
  M        <- 3 
  nvec     <- rep(500, 3)
  n        <- sum(nvec)
  
  lambda   <- 0.4              # peer effects
  beta     <- c(2, -1.9, 0.8)  # own effects
  gamma    <- c(1.5, -1.2)     # contextual effects
  sigma    <- 1.5              # standard deviation of epsilon
  theta    <- c(lambda, beta, gamma, sigma)
  
  data     <- data.frame(x1 = rnorm(n, 1, 1), x2 =  rexp(n, 0.4))
  
  # Network
  betanet  <- c(-2.8, -1.5)    # beta
  Glist    <- list()           # adjacency matrix row normalized
  Network  <- list()           # adjacency matrix row non-normalized
  dX       <- matrix(0, 0, 2)  # observed dyad-specific variables
  mu       <- list()           # unobserved individual-level attribute
  uu       <- runif(M, -1, 1)  # mean of uu in each sub-network
  sigma2u  <- runif(M, 0.5, 4) # variance of uu in each sub-network
  
  for (m in 1:M) {
    nm           <- nvec[m]
    mum          <- rnorm(nm, uu[m], sqrt(sigma2u[m]))
    Z1           <- matrix(0, nm, nm)  
    Z2           <- matrix(0, nm, nm)
    
    for (i in 1:nm) {
      for (j in 1:nm) {
        Z1[i, j] <- abs(data$x1[i] - data$x1[j])
        Z2[i, j] <- abs(data$x2[i] - data$x2[j])
      }
    }
    
    Gm           <- 1*((Z1*betanet[1] + Z2*betanet[2] +
                          kronecker(mum, t(mum), "+") + rlogis(nm^2)) > 0)
    diag(Gm)     <- 0
    
    diag(Z1)     <- NA
    diag(Z2)     <- NA
    Z1           <- Z1[!is.na(Z1)]
    Z2           <- Z2[!is.na(Z2)]
    
    dX           <- rbind(dX, cbind(Z1, Z2))
    
    Network[[m]] <- Gm
    rs           <- rowSums(Gm); rs[rs == 0] <- 1
    Gm           <- Gm/rs
    Glist[[m]]   <- Gm
    mu[[m]]      <- mum
  }
  mu             <- unlist(mu)
  
  # Data
  tmu      <- (mu - rep(uu, nvec))/sqrt(rep(sigma2u, nvec))
  data$tmu <- tmu
  rho      <- 0.24
  rhobar   <- 0.18
  thetanet <- c(lambda, beta, sigma*rho, gamma, sigma*rhobar,
                sigma*(1 - rho - rhobar))
  ytmp     <- simCDnet(formula = ~ x1 + x2 + tmu | x1 + x2 + tmu,
                       Glist = Glist, theta = thetanet, data = data)
  y        <- ytmp$y
  data$y   <- y
  
  # Count data model
  CDexo <- CDnetNPL(formula = y ~ x1 + x2, contextual = TRUE, Glist = Glist,
                    optimizer = "optim", data = data, npl.ctr = list(print = FALSE))
  # Dyadic linking model
  net   <- netformation(network =  Network, formula = ~ dX, fixed.effects = TRUE,
                        mcmc.ctr = list(burnin = 1000, iteration = 2000), print = FALSE)
  # Endogeneity
  t           <- which.max(net$posterior$log.density)
  muest       <- net$posterior$mu[t,]
  uuest       <- net$posterior$uu[t,]
  sigma2uest  <- net$posterior$sigmamu2[t,]
  tmuest      <- (muest - rep(uuest, nvec))/sqrt(rep(sigma2uest, nvec))
  data$tmuest <- tmuest
  CDendo      <- CDnetNPL(formula = y ~ x1 + x2 + tmuest, contextual = TRUE, 
                          Glist = Glist, optimizer = "optim", data = data, 
                          npl.ctr = list(print = FALSE))
  c(CDexo$estimate, CDendo$estimate)
}

## ----mcendoout, echo = TRUE, eval = FALSE, message = FALSE--------------------
#  n.cores               <- 5
#  replic                <- 1000 #Number of replications
#  out.mcendo            <- mclapply(1:replic, fMCendo, mc.cores = n.cores)
#  out.mcendo            <- apply(t(do.call(cbind, out.mcendo)), 2, mean)
#  out.endo              <- cbind(thetanet, c(out.mcendo[c(1:4, NA, 5:6, NA, 7)]),
#                                 out.mcendo[8:16])
#  rownames(out.endo)    <- names(out.mcendo[8:16])
#  colnames(out.endo)    <- c("TrueValue", "Endo-Not-Controlled", "Endo-Controlled")
#  print(out.endo)

## ----mcprintendo, echo = FALSE, eval = TRUE-----------------------------------
load("out.mcendo.Rdata")
print(out.endo)

