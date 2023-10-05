library(Rcpp)
library(tidyverse)
sourceCpp("mosumCwald.cpp")

MOSUM.cp.candidate <- function(y, time, G, eta = .3, alpha = .05){
  # y : vector of observations
  # G_vec : Vector of bandwidths to be used 
  # eta : a numeric value in (0,1/2) for the minimal size of exceeding environments

  cp <- vector()
  n <- length(y)
  
  # calculate threshold 
  logH <- 0.7284
  b0_G <- 2*log(n/G) + log(log(n/G)) + logH
  a_G <- sqrt(2*log(n/G))
  q_G <- (b0_G - log(-log(1 - alpha)/2))/a_G      
  
  # Wald-type statistics
  TG <- mosumCwald(y, time, G=G)
  
  # Collect intervals that the test statistics exceeds threshold
  reduced.time.index <- c(G:(n-G))
  stat.over.thres <- (TG >= q_G)
  poo <- rle(stat.over.thres)
  criterion <- (poo$lengths >= (eta*G))
  index.over.thres <- cumsum(poo$lengths)
  search.cand <- which(poo$values & criterion)
  cp_G <- TG_G <- vector()
  if(length(search.cand)>0){
    for(j in 1:length(search.cand)){
      i <- search.cand[j]
      if(i!=1){
        search.cand.loop <- c((index.over.thres[i-1]+1):(index.over.thres[i]))
        cp_G[j] <- search.cand.loop[which.max(TG[search.cand.loop])]
        TG_G[j] <- max(TG[search.cand.loop])
      }else{
        search.cand.loop <- c(1:(index.over.thres[i]))
        cp_G[j] <- search.cand.loop[which.max(TG[search.cand.loop])]
        TG_G[j] <- max(TG[search.cand.loop])
      }
      TG_G <- TG_G[!is.na(cp_G)]
      cp_G <- cp_G[!is.na(cp_G)]
    }
  }

  cp_G.index <- reduced.time.index[cp_G]
  
  return(data.frame(cp = cp_G.index, TG = TG_G))
}

MOSUM_linear <- function(y, time=c(1:length(y)), G_vec, eta=.3, theta=.8, alpha=.05, 
                         sort.G = c("bottomup", "aic", "bic", "cvc")[3]){
  # y : vector of observations
  # G_vec : Vector of bandwidths to be used 
  # eta : a numeric value in (0,1/2) for the minimal size of exceeding environments
  # theta : a numeric value in (0,1] for the minimal size of accepting new candidates in the merging procedure
  # alpha : significance level
  # sort.G : how to sort bandwidths? 
  
  cp <- vector()
  n <- length(y)
  G_vec.sort <- sort(G_vec, decreasing=FALSE)
  cp_cand_df <- data.frame(G=NULL, cp=NULL, wald=NULL)
  
  for(lag.idx in 1:length(G_vec.sort)){
    lag <- G_vec.sort[lag.idx]
    
    # find candidates
    cp.candidate.obj <- MOSUM.cp.candidate(y, time, G = lag, eta = eta, alpha = alpha)
    cp_G.index <- cp.candidate.obj$cp
    TG_G <- cp.candidate.obj$TG
    
    # Obtain residual sum of squares
    if(length(cp_G.index)>0){
      rss <- mosumCwald_rss(y, cp_index=cp_G.index)
    }else{
      rss <- slr_rss(y)
    }
    
    # data frame with bandwidths and change points candidates
    if(sort.G == "cvc"){
      odd <- seq(from = 1, to = n, by = 2)
      even <- seq(from = 2, to = n, by = 2)
      cp_odd  <- MOSUM.cp.candidate(y[odd], time[odd], G = floor(lag/2), eta = eta, crit = crit, alpha = alpha)$cp
      cp_even <- MOSUM.cp.candidate(y[even], time[even], floor(lag/2), eta = eta, crit = crit, alpha = alpha)$cp
      cvc <- cv_criterion(data = y, train = odd, test = even, cp = cp_odd) + 
        cv_criterion(data = y, train = even, test = odd, cp = cp_even)
      cp_cand_df <- rbind(cp_cand_df, data.frame(G=rep(lag, length(cp_G.index)), cp=cp_G.index, 
                                                 cp.num=rep(length(cp_G.index), length(cp_G.index)),
                                                 wald=TG_G, rss=rep(rss, length(cp_G.index)), 
                                                 cvc = rep(cvc, length(cp_G.index))
      ))
    }else{
      cp_cand_df <- rbind(cp_cand_df, data.frame(G=rep(lag, length(cp_G.index)), cp=cp_G.index, 
                                                 cp.num=rep(length(cp_G.index), length(cp_G.index)),
                                                 wald=TG_G, rss=rep(rss, length(cp_G.index))
      ))
    }
    
  }
  
  # obtain AIC and BIC 
  aic <- n*log(cp_cand_df$rss/n) + 2*2*(cp_cand_df$cp.num+1)
  bic <- n*log(cp_cand_df$rss/n) + 2*log(n)*(cp_cand_df$cp.num+1)
  cp_cand_df$aic <- aic
  cp_cand_df$bic <- bic

  # Merging procedure
  if(sort.G == "bottomup"){
    for(i in seq(from=1, to=dim(cp_cand_df)[1])){
      cp_cand <- cp_cand_df$cp[i]
      G <- cp_cand_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else if(sort.G == "bic"){
    search_ordered_df <- cp_cand_df %>% group_by(G) %>% arrange(bic, desc(wald)) %>% ungroup()
    for(i in 1:dim(search_ordered_df)[1]){           
      cp_cand <- search_ordered_df$cp[i]
      G <- search_ordered_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else if(sort.G == "aic"){
    search_ordered_df <- cp_cand_df %>% group_by(G) %>% arrange(aic, desc(wald)) %>% ungroup()
    for(i in 1:dim(search_ordered_df)[1]){           
      cp_cand <- search_ordered_df$cp[i]
      G <- search_ordered_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else if(sort.G == "cvc"){
    search_ordered_df <- cp_cand_df %>% group_by(G) %>% arrange(cvc, desc(wald)) %>% ungroup()
    for(i in 1:dim(search_ordered_df)[1]){           
      cp_cand <- search_ordered_df$cp[i]
      G <- search_ordered_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else{
    print("Error: Wrong input of sorting argument")
    break
  }
  return(sort(cp))
}


### Dependent error 

# LRV estimator by Chan (2022)
# This code comes from the author's homepage 
dlrv_simple = function(x){
  n = length(x)
  l = ceiling(2*n^(1/5))
  h = 2*l
  m = 3
  d = c(0.1942, 0.2809, 0.3832, -0.8582)
  D = rep(0,n-m*h)
  for(j in 0:m) D = D + d[j+1] * x[(m*h+1-j*h):(n-j*h)]
  acvf = acf(D,type="cov",plot=FALSE ,lag.max=l)$acf
  acvf[1] + 2*sum(acvf[-1]*(1-((1:l)/l)^2))
}

# Candidate using LRV estimator

MOSUM.cp.candidate.ar <- function(y, time, G, eta = .3, alpha = .05){
  # y : vector of observations
  # G_vec : Vector of bandwidths to be used 
  # eta : a numeric value in (0,1/2) for the minimal size of exceeding environments
  
  cp <- vector()
  n <- length(y)

  # calculate threshold 
  logH <- 0.7284
  b0_G <- 2*log(n/G) + log(log(n/G)) + logH
  a_G <- sqrt(2*log(n/G))
  q_G <- (b0_G - log(-log(1 - alpha)/2))/a_G    
  
  # Wald-type statistics
  
  TG <- mosumCwald_nonsigma(y, time, G=G)

  lrv_vec <- vector(length = length(y) - G + 1)
  for(i in 1:length(lrv_vec)){
    lrv_vec[i] <- dlrv_simple(y[i:(i+G-1)])
  }
  lrv <- (lrv_vec[(G+1):length(lrv_vec)] + lrv_vec[1:(length(lrv_vec)-G)])/2
  TG <- TG / sqrt(max(lrv, 1e-04))
  
  # Collect intervals that the test statistics exceeds threshold
  reduced.time.index <- c(G:(n-G))
  stat.over.thres <- (TG >= q_G)
  poo <- rle(stat.over.thres)
  criterion <- (poo$lengths >= (eta*G))
  index.over.thres <- cumsum(poo$lengths)
  search.cand <- which(poo$values & criterion)
  cp_G <- TG_G <- vector()
  if(length(search.cand)>0){
    for(j in 1:length(search.cand)){
      i <- search.cand[j]
      if(i!=1){
        search.cand.loop <- c((index.over.thres[i-1]+1):(index.over.thres[i]))
        cp_G[j] <- search.cand.loop[which.max(TG[search.cand.loop])]
        TG_G[j] <- max(TG[search.cand.loop])
      }else{
        search.cand.loop <- c(1:(index.over.thres[i]))
        cp_G[j] <- search.cand.loop[which.max(TG[search.cand.loop])]
        TG_G[j] <- max(TG[search.cand.loop])
      }
      TG_G <- TG_G[!is.na(cp_G)]
      cp_G <- cp_G[!is.na(cp_G)]
    }
  }
  
  cp_G.index <- reduced.time.index[cp_G]
  
  return(data.frame(cp = cp_G.index, TG = TG_G))
}

MOSUM_linear_ar <- function(y, time, G_vec, eta=.3, theta=.8, 
                                       sort.G = c("bottomup", "aic", "bic", "cvc")[3]){
  # y : vector of observations
  # G_vec : Vector of bandwidths to be used 
  # eta : a numeric value in (0,1/2) for the minimal size of exceeding environments
  # theta : a numeric value in (0,1] for the minimal size of accepting new candidates in the merging procedure
  # alpha : significance level
  # sort.G : how to sort bandwidths? 
  
  cp <- vector()
  n <- length(y)
  G_vec.sort <- sort(G_vec, decreasing=FALSE)
  cp_cand_df <- data.frame(G=NULL, cp=NULL, wald=NULL)
  for(lag.idx in 1:length(G_vec.sort)){
    lag <- G_vec.sort[lag.idx]

    # find candidates
    cp.candidate.obj <- MOSUM.cp.candidate.ar(y, time=c(1:length(y)), G = lag, eta = eta)
    cp_G.index <- cp.candidate.obj$cp
    TG_G <- cp.candidate.obj$TG
    
    # Obtain residual sum of squares
    if(length(cp_G.index)>0){
      rss <- mosumCwald_rss(y, cp_index=cp_G.index)
    }else{
      rss <- slr_rss(y)
    }

    # data frame with bandwidths and change points candidates
    if(sort.G == "cvc"){
      odd <- seq(from = 1, to = n, by = 2)
      even <- seq(from = 2, to = n, by = 2)
      cp_odd  <- MOSUM.cp.candidate.ar(y[odd], time[odd], G = floor(lag/2), eta = eta, crit = crit)$cp
      cp_even <- MOSUM.cp.candidate.ar(y[even], time[even], floor(lag/2), eta = eta, crit = crit)$cp
      cvc <- cv_criterion(data = y, train = odd, test = even, cp = cp_odd) + 
        cv_criterion(data = y, train = even, test = odd, cp = cp_even)
      cp_cand_df <- rbind(cp_cand_df, data.frame(G=rep(lag, length(cp_G.index)), cp=cp_G.index, 
                                                 cp.num=rep(length(cp_G.index), length(cp_G.index)),
                                                 wald=TG_G, rss=rep(rss, length(cp_G.index)), 
                                                 cvc = rep(cvc, length(cp_G.index))
      ))
    }else{
      cp_cand_df <- rbind(cp_cand_df, data.frame(G=rep(lag, length(cp_G.index)), cp=cp_G.index, 
                                                 cp.num=rep(length(cp_G.index), length(cp_G.index)),
                                                 wald=TG_G, rss=rep(rss, length(cp_G.index))
      ))
    }
  }

  # Obtain AIC and bic
  aic <- n*log(cp_cand_df$rss/n) + 2*2*(cp_cand_df$cp.num+1)
  bic <- n*log(cp_cand_df$rss/n) + 2*log(n)*(cp_cand_df$cp.num+1)
  cp_cand_df$aic <- aic
  cp_cand_df$bic <- bic

  # Merging procedure
  if(sort.G == "bottomup"){
    for(i in seq(from=1, to=dim(cp_cand_df)[1])){
      cp_cand <- cp_cand_df$cp[i]
      G <- cp_cand_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else if(sort.G == "bic"){
    search_ordered_df <- cp_cand_df %>% group_by(G) %>% arrange(bic, desc(wald)) %>% ungroup()
    for(i in 1:dim(search_ordered_df)[1]){           
      cp_cand <- search_ordered_df$cp[i]
      G <- search_ordered_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else if(sort.G == "aic"){
    search_ordered_df <- cp_cand_df %>% group_by(G) %>% arrange(aic, desc(wald)) %>% ungroup()
    for(i in 1:dim(search_ordered_df)[1]){           
      cp_cand <- search_ordered_df$cp[i]
      G <- search_ordered_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else if(sort.G == "cvc"){
    search_ordered_df <- cp_cand_df %>% group_by(G) %>% arrange(cvc, desc(wald)) %>% ungroup()
    for(i in 1:dim(search_ordered_df)[1]){           
      cp_cand <- search_ordered_df$cp[i]
      G <- search_ordered_df$G[i]
      if(length(cp_cand)>0){
        if(length(cp)==0){
          cp <- c(cp, cp_cand)
        }else if(min(abs(cp_cand-cp))>=theta*G){
          cp <- c(cp, cp_cand)
        }
      }
    }
  }else{
    print("Error: Wrong input of sorting argument")
    break
  }
  return(sort(cp))
}
