source("MOSUM_linear.R")

## Example (True change points at 1000, 2000, 2500)
beta <- c(-1, -1, -2.5, 2.5)+rnorm(4, sd=0.2)
underlying <- c(seq(beta[1]*(0.01-10)+10, 10, length=1000), beta[2]*seq(from=0.01, to=10, length=1000), 
                10*(1+beta[2])+beta[3]*seq(from=0.01, to=5, length=500), 
                10*(1+beta[2])+5*beta[3]+beta[4]*seq(from=0.01, to=10, length=1000))
y<-underlying+rnorm(length(underlying),0,sd=1)

G_vec <- c(50, 100, 150, 250, 400, 650)
MOSUM_linear(y, G_vec=G_vec, sort.G="bic")


## Simulation 1: (M1) with (E1)
set.seed(191009)
cp.true <- c(10, 20, 25)
simuldata <- vector('list', length=6)
sd.vec <- c(0.5, 1, 1.5, 2)
for(snr in 1:4){
  sd.noise <- sd.vec[snr]
  simuldata[[snr]] <- vector('list', length=1000)
  for(iter in 1:1000){
    beta <- c(-1, -1, -2.5, 2.5)+rnorm(4, sd=0.2)
    underlying <- c(seq(beta[1]*(0.01-10)+10, 10, length=1000), beta[2]*seq(from=0.01, to=10, length=1000), 10*(1+beta[2])+beta[3]*seq(from=0.01, to=5, length=500), 10*(1+beta[2])+5*beta[3]+beta[4]*seq(from=0.01, to=10, length=1000))
    y<-underlying+rnorm(length(underlying),0,sd=sd.noise)
    time <- (1:length(y))/100
    simuldata[[snr]][[iter]]$time <- time
    simuldata[[snr]][[iter]]$value <- y
    simuldata[[snr]][[iter]]$underlying <- underlying
  }
}

cp.estimated <- vector('list', length=7)
names(cp.estimated) <- c("sigma=0.5", "sigma=1", "sigma=1.5", "sigma=2", "params", "time", "description")
for(snr in 1:4){
  cp.estimated[[snr]] <- vector('list', length=1000)
}

G_vec <- c(50, 100, 150, 250, 400, 650)
eta.simul <- .3
theta.simul <- 0.8

for(snr in 1:4){
  for(iter in 1:1000){
    # 
    # 
    
    y <- simuldata[[snr]][[iter]]$value
    time <- simuldata[[snr]][[iter]]$time
    
    poo.cp <- MOSUM_linear(y, time, G_vec, eta=eta.simul, theta=theta.simul, alpha=.05, sort.G="bic")
    cp.estimated[[snr]][[iter]] <- time[poo.cp]
  }
}

cp.estimated$params <- c(eta.simul, theta.simul)
names(cp.estimated$params) <- c("eta", "theta")
cp.estimated$time <- mosum.time
cp.estimated$description <- "Simulation 1, 1000 iterations, 6 bandwidths (G=50, 100, 150, 250, 400, 650), all 4 sigmas, BIC-Up Merging"

### Performance evaluation
maxscore.L1 <- function(x, cp.true){
  if(length(x)==0){
    return(Inf)
  }else{
    return.val <- vector(length=length(cp.true))
    for(j in 1:length(cp.true)){
      return.val[j] <- min(abs(x - cp.true[j]))
    }
    return(max(return.val))
  }
}
count.L1 <- function(x, cp){
  return(abs(length(x)-length(cp)))
}
maxscore.inv <- function(x, cp.true){
  if(length(x)==0){
    return(-Inf)
  }else{
    return(maxscore.L1(cp.true, x))
  }
}

cp.true <- c(10, 20, 25)
result.ls.mosum <- vector('list', length=4)
names(result.ls.mosum) <- c("sigma=0.5", "sigma=1", "sigma=1.5", "sigma=2")
for(snr in 1:4){
  result.ls.mosum[[snr]] <- vector('list', length=3)
  names(result.ls.mosum[[snr]]) <- c("MAXscore", "MAXscore2", "COUNTscore")
  result.ls.mosum[[snr]]$MAXscore <- unlist(lapply(cp.estimated[[snr]], maxscore.L1, cp.true=cp.true))
  result.ls.mosum[[snr]]$MAXscore2 <- unlist(lapply(cp.estimated[[snr]], maxscore.inv, cp.true=cp.true))
  result.ls.mosum[[snr]]$COUNTscore <- unlist(lapply(cp.estimated[[snr]], count.L1, cp=cp.true))
}

maxscore_mosum_simul1 <- array(dim=c(4, 1000))
maxscore2_mosum_simul1 <- array(dim=c(4, 1000))
count_mosum_simul1 <- array(dim=c(4, 1000))

for(snr in 1:4){
  maxscore_mosum_simul1[snr,] <- result.ls.mosum[[snr]]$MAXscore
  maxscore2_mosum_simul1[snr,] <- result.ls.mosum[[snr]]$MAXscore2
  count_mosum_simul1[snr,] <- result.ls.mosum[[snr]]$COUNTscore
}
hausdorff_mosum_simul1 <- pmax(maxscore_mosum_simul1, maxscore2_mosum_simul1)

apply(count_mosum_simul1, 1, mean)
apply(maxscore_mosum_simul1, 1, mean)
apply(maxscore2_mosum_simul1, 1, mean)
apply(hausdorff_mosum_simul1, 1, mean)


## Simulation 2: (M2) with (E2)
set.seed(220518)
cp.true <- c(10, 20, 25)
simuldata2 <- vector('list', length=6)
sd.vec <- c(0.5, 1, 1.5, 2)
for(snr in 1:4){
  sd.noise <- sd.vec[snr]
  simuldata2[[snr]] <- vector('list', length=1000)
  for(iter in 1:1000){
    beta <- c(-1, 1, -2.5, 2.5)+rnorm(4, sd=0.2)
    underlying <- c(seq(beta[1]*(0.01-10), 0, length=1000), beta[2]*seq(from=0.01, to=10, length=1000), 10*beta[2]+beta[3]*seq(from=0.01, to=5, length=500), 10*beta[2]+5*beta[3]+beta[4]*seq(from=0.01, to=10, length=1000))
    y<-underlying+rnorm(length(underlying),0,sd=sd.noise)
    time <- (1:length(y))/100
    simuldata2[[snr]][[iter]]$time <- time
    simuldata2[[snr]][[iter]]$value <- y
    simuldata2[[snr]][[iter]]$underlying <- underlying
  }
}

## Simulation 3: (M1) under (E4)
set.seed(230602)
cp.true <- c(10, 20, 25)
rho <- 0.3
# rho <- 0.7
simuldata <- vector('list', length=4)
snr.vec <- c(0.5, 1, 1.5, 2) * sqrt(1-rho^2) # this is sigma_delta. 
for(snr in 1:4){
  sd.noise <- snr.vec[snr] # this is sigma_delta. 
  simuldata[[snr]] <- vector('list', length=1000)
  for(iter in 1:1000){
    
    beta <- c(-1, -1, -2.5, 2.5)+rnorm(4, sd=0.2)
    underlying <- c(seq(beta[1]*(0.01-10)+10, 10, length=1000), beta[2]*seq(from=0.01, to=10, length=1000), 10*(1+beta[2])+beta[3]*seq(from=0.01, to=5, length=500), 10*(1+beta[2])+5*beta[3]+beta[4]*seq(from=0.01, to=10, length=1000))
    # len_eps <- length(underlying) + 500 
    # epsilon <- vector(length = len_eps)
    # epsilon[1] <- rnorm(1, sd = sd.noise)
    # for(i in 2:len_eps){
    #   epsilon[i] <- rho * epsilon[i-1] + rnorm(1, sd = sd.noise)
    # }
    # epsilon <- epsilon[501:len_eps]
    epsilon <- arima.sim(n = 3500, list(ar = rho), sd = sd.noise)
    y<-underlying+epsilon
    time <- (1:length(y))/100
    simuldata[[snr]][[iter]]$time <- time
    simuldata[[snr]][[iter]]$value <- y
    simuldata[[snr]][[iter]]$underlying <- underlying
  }
}

y <- simuldata[[3]][[1]]$value
G_vec <- c(50, 100, 150, 250, 400, 650)
MOSUM_linear(y, G_vec=G_vec, sort.G="bic") # MOSUM
MOSUM_linear_ar(y, G_vec=G_vec, sort.G="bic") # MOSUMdlrv
