library(Rcpp)
sourceCpp("mosumCwald.cpp")

MOSUM_linear <- function(y, time=c(1:length(y)), G_vec, eta=.3, theta=.8, alpha=.05, sort.G = c("default", "aic", "bic")[1]){
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
    
    # calculate threshold 
    logH <- 0.7284
    b0_G <- 2*log(n/lag) + log(log(n/lag)) + logH
    a_G <- sqrt(2*log(n/lag))
    q_G <- (b0_G - log(-log(1-alpha)/2))/a_G  
    
    # Wald-type statistics
    TG <- mosumCwald(y, time, G=lag)
    
    # Collect intervals that the test statistics exceeds threshold
    reduced.time.index <- c(lag:(n-lag))
    stat.over.thres <- (TG >= q_G)
    poo <- rle(stat.over.thres)
    criterion <- (poo$lengths >= (eta*lag))
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
    
    # Obtain residual sum of squares
    cp_G.index <- reduced.time.index[cp_G]
    if(length(cp_G.index)>0){
      rss <- mosumCwald_rss(y, cp_index=cp_G.index)
    }else{
      rss <- slr_rss(y)
    }
    
    # data frame with bandwidths and change points candidates
    cp_cand_df <- rbind(cp_cand_df, data.frame(G=rep(lag, length(cp_G.index)), cp=cp_G.index, 
                                               cp.num=rep(length(cp_G.index), length(cp_G.index)),
                                               wald=TG_G, rss=rep(rss, length(cp_G.index))))
  }
  
  # obtain AIC and BIC 
  aic <- n*log(cp_cand_df$rss/n) + 2*2*(cp_cand_df$cp.num+1)
  bic <- n*log(cp_cand_df$rss/n) + 2*log(n)*(cp_cand_df$cp.num+1)
  cp_cand_df$aic <- aic
  cp_cand_df$bic <- bic

  # Merging procedure
  if(sort.G == "default"){
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
    searching_order <- order(cp_cand_df$bic, decreasing=FALSE)
    for(i in searching_order){           
      cp_cand <- cp_cand_df$cp[i]
      G <- cp_cand_df$G[i]
      if(length(cp)==0){
        cp <- c(cp, cp_cand)
      }else if(min(abs(cp_cand-cp))>=theta*G){
        cp <- c(cp, cp_cand)
      }
    }
  }else if(sort.G == "aic"){
    searching_order <- order(cp_cand_df$aic, decreasing=FALSE)
    for(i in searching_order){           
      cp_cand <- cp_cand_df$cp[i]
      G <- cp_cand_df$G[i]
      if(length(cp)==0){
        cp <- c(cp, cp_cand)
      }else if(min(abs(cp_cand-cp))>=theta*G){
        cp <- c(cp, cp_cand)
      }
    }
  }else{
    print("Error: Wrong input of sorting argument")
    break
  }
  return(sort(cp))
}
