# Change Points Detection based on Moving Sum Method under Piecewise Linearity

R function `MOSUM_linear` conducts the moving sum based change points detection proposed by 
> Kim, J., Oh, H. S., & Cho, H. (2022+). <a href = "https://arxiv.org/abs/2208.04900">Moving sum procedure for change point detection under piecewise linearity</a>. 

To implement `MOSUM_linear` function, you first have to source `MOSUM_linear.R` file in R. It includes the procedure loading functions in `mosumCwald.cpp`. For this `Rcpp` package should be preliminarily installed.
```{r}
source("MOSUM_linear.R")
```

`MOSUM_linear` function requires following input arguments. 

- `y`: a data vector to be segmented.
- `time`: a time vector. Default is `c(1:length(y))`.
- `G_vec`: a vector of bandwidths.
- `eta`: a value of $\eta$. For details see page 5 of the paper. Default value is 0.3.
- `theta`: a value of $\theta$. For details see Section 3.2 of the paper. Default value is 0.8.
- `alpha`: significance level. Default is 0.05.
- `sort.G`: a merge criterion. `bottomup`, `aic`, `bic`, `cvc` are available. Default is `bic`. 

Following code illustrates a toy example. 

```{r}
beta <- c(-1, -1, -2.5, 2.5)+rnorm(4, sd=0.2)
underlying <- c(seq(beta[1]*(0.01-10)+10, 10, length=1000), beta[2]*seq(from=0.01, to=10, length=1000), 
                10*(1+beta[2])+beta[3]*seq(from=0.01, to=5, length=500), 
                10*(1+beta[2])+5*beta[3]+beta[4]*seq(from=0.01, to=10, length=1000))
y<-underlying+rnorm(length(underlying),0,sd=1)

G_vec <- c(50, 100, 150, 250, 400, 650)
MOSUM_linear(y, G_vec=G_vec, sort.G="bic")
```

`simulations.R` file conducts a part of simulation demonstrated in the manuscript. 
