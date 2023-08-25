// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <random>
#include <math.h>
#include <cmath>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mosumCwald(NumericVector data, NumericVector time, int G) {
    int n=data.size();
    NumericVector out(n-2*G+1);
    double dt = time(n-1) / n;
    double G_flt = double(G);
    double G_half = (G_flt+1.)/2.;
    
    arma::vec diff_data(n-G);
    
    for(int i=0; i<n-G; i++){
        diff_data(i) = data(G+i) - data(i);
    }
    
    double denom = 0;
    double beta1plus = 0;
    double beta1minus = 0;
    double beta0diff = 0;
    
    /* Initializa */
    for(int i=0; i<G; i++){
        double i_flt = double(i+1);
        double w = i_flt - G_half;
        denom += pow(w, 2.0);
        beta1plus += w*data(G+i);
        beta1minus += w*data(i);
        beta0diff += diff_data(i);
    }
    
    denom *= dt;
    
    beta1plus /= denom;
    beta1minus /= denom;
    beta0diff = beta0diff/G_flt - (G_half*beta1plus + (G_half-1.)*beta1minus)*dt;
    
    double sum_plus = 0; double sum_minus = 0;
    double sumsq_plus = 0; double sumsq_minus = 0;
    for(int i=0; i<G; i++){
        sum_plus += data(G+i);
        sum_minus += data(i);
        sumsq_plus += pow(data(G+i), 2.);
        sumsq_minus += pow(data(i), 2.);
    }
    
    double SSTplus = sumsq_plus - pow(sum_plus, 2.)/G_flt;
    double SSTminus = sumsq_minus - pow(sum_minus, 2.)/G_flt;
    double sigmasq_plus = (SSTplus - pow(beta1plus, 2.)*denom*dt)/(G_flt-2.);
    double sigmasq_minus = (SSTminus - pow(beta1minus, 2.)*denom*dt)/(G_flt-2.);
    
    double sigma_sq = (sigmasq_plus + sigmasq_minus)/2.;
    
    out(0) = sqrt(G_flt/(8.*sigma_sq)*(pow(beta0diff, 2.) + pow(G_flt*(beta1plus - beta1minus)*dt, 2.)/3.));
    
    /* update */
    double Delta_b1plus = 0;
    double Delta_b1minus = 0;
    double Delta_b0diff = 0;
    for(int k=G-1; k<n-G-1; k++){
        Delta_b1plus = (G_flt*data(k+G+1) - sum_plus - G_half*diff_data(k+1))/denom;
        Delta_b1minus = (G_flt*data(k+1) - sum_minus - G_half*diff_data(k-G+1))/denom;
        Delta_b0diff = (diff_data(k+1) - diff_data(k-G+1))/G_flt - (G_half*Delta_b1plus + (G_half-1.)*Delta_b1minus)*dt;
        beta1plus += Delta_b1plus;
        beta1minus += Delta_b1minus;
        beta0diff += Delta_b0diff;
        sum_plus += data(k+G+1) - data(k+1);
        sum_minus += data(k+1) - data(k-G+1);
        sumsq_plus += pow(data(k+G+1), 2.) - pow(data(k+1), 2.);
        sumsq_minus += pow(data(k+1), 2.) - pow(data(k-G+1), 2.);
        SSTplus = sumsq_plus - pow(sum_plus, 2.)/G_flt;
        SSTminus = sumsq_minus - pow(sum_minus, 2.)/G_flt;
        sigmasq_plus = (SSTplus - pow(beta1plus, 2.)*denom*dt)/(G_flt-2.);
        sigmasq_minus = (SSTminus - pow(beta1minus, 2.)*denom*dt)/(G_flt-2.);
        sigma_sq = (sigmasq_plus + sigmasq_minus)/2.;
        out(k-G+2) = sqrt(G_flt/8.*(pow(beta0diff, 2.) + pow(G_flt*(beta1plus - beta1minus)*dt, 2.)/3.)/sigma_sq);
    }
    
    return out;
    
}


// [[Rcpp::export]]
NumericVector mosumCwald_var(NumericVector data, NumericVector time, int G) {
    int n=data.size();
    NumericVector out(n-2*G+1);
    double dt = time(n-1) / n;
    double G_flt = double(G);
    double G_half = (G_flt+1.)/2.;
    
    arma::vec diff_data(n-G);
    
    for(int i=0; i<n-G; i++){
        diff_data(i) = data(G+i) - data(i);
    }
    
    double denom = 0;
    double beta1plus = 0;
    double beta1minus = 0;
    double beta0diff = 0;
    
    /* Initializa */
    for(int i=0; i<G; i++){
        double i_flt = double(i+1);
        double w = i_flt - G_half;
        denom += pow(w, 2.0);
        beta1plus += w*data(G+i);
        beta1minus += w*data(i);
        beta0diff += diff_data(i);
    }
    
    denom *= dt;
    
    beta1plus /= denom;
    beta1minus /= denom;
    beta0diff = beta0diff/G_flt - (G_half*beta1plus + (G_half-1.)*beta1minus)*dt;
    
    double sum_plus = 0; double sum_minus = 0;
    double sumsq_plus = 0; double sumsq_minus = 0;
    for(int i=0; i<G; i++){
        sum_plus += data(G+i);
        sum_minus += data(i);
        sumsq_plus += pow(data(G+i), 2.);
        sumsq_minus += pow(data(i), 2.);
    }
    
    double SSTplus = sumsq_plus - pow(sum_plus, 2.)/G_flt;
    double SSTminus = sumsq_minus - pow(sum_minus, 2.)/G_flt;
    double sigmasq_plus = (SSTplus - pow(beta1plus, 2.)*denom*dt)/(G_flt-2.);
    double sigmasq_minus = (SSTminus - pow(beta1minus, 2.)*denom*dt)/(G_flt-2.);

    double sigma_sq = (sigmasq_plus + sigmasq_minus)/2.;
    
    out(0) = sigma_sq;
    
    /* update */
    double Delta_b1plus = 0;
    double Delta_b1minus = 0;
    double Delta_b0diff = 0;
    for(int k=G-1; k<n-G-1; k++){
        Delta_b1plus = (G_flt*data(k+G+1) - sum_plus - G_half*diff_data(k+1))/denom;
        Delta_b1minus = (G_flt*data(k+1) - sum_minus - G_half*diff_data(k-G+1))/denom;
        Delta_b0diff = (diff_data(k+1) - diff_data(k-G+1))/G_flt - (G_half*Delta_b1plus + (G_half-1.)*Delta_b1minus)*dt;
        beta1plus += Delta_b1plus;
        beta1minus += Delta_b1minus;
        beta0diff += Delta_b0diff;
        sum_plus += data(k+G+1) - data(k+1);
        sum_minus += data(k+1) - data(k-G+1);
        sumsq_plus += pow(data(k+G+1), 2.) - pow(data(k+1), 2.);
        sumsq_minus += pow(data(k+1), 2.) - pow(data(k-G+1), 2.);
        SSTplus = sumsq_plus - pow(sum_plus, 2.)/G_flt;
        SSTminus = sumsq_minus - pow(sum_minus, 2.)/G_flt;
        sigmasq_plus = (SSTplus - pow(beta1plus, 2.)*denom*dt)/(G_flt-2.);
        sigmasq_minus = (SSTminus - pow(beta1minus, 2.)*denom*dt)/(G_flt-2.);
        sigma_sq = (sigmasq_plus + sigmasq_minus)/2.;
        out(k-G+2) = sigma_sq;
    }

    return out;

}

// [[Rcpp::export]]
NumericVector mosumCwald_fit(NumericVector data, NumericVector cp_index) {
    int n=data.size();
    int J=cp_index.size();
    NumericVector out(n);
    
    double Sxx=0;
    double Sxy=0;
    double xbar=0;
    double ybar=0;
        
    for(int i=0; i<cp_index(0); i++){
        xbar += i;
        ybar += data(i);
        Sxx += double(pow(i, 2.));
        Sxy += double(i)*data(i);
    }
    xbar = double(xbar) / cp_index(0);
    ybar = ybar / cp_index(0);
    Sxx = Sxx - cp_index(0)*pow(double(xbar), 2.);
    Sxy = Sxy - cp_index(0)*double(xbar)*ybar;
    
    
    for(int i=0; i<cp_index(0); i++){
        out(i) = ybar + Sxy/Sxx * (double(i)-xbar);
    }
    
    
    for(int j=0; j<J-1; j++){
        Sxx=0;
        Sxy=0;
        xbar=0;
        ybar=0;
        
        for(int i=cp_index(j); i<cp_index(j+1); i++){
            xbar += i;
            ybar += data(i);
            Sxx += double(pow(i, 2.));
            Sxy += double(i)*data(i);
        }
        xbar = double(xbar) / (cp_index(j+1) - cp_index(j));
        ybar = ybar / (cp_index(j+1) - cp_index(j));
        Sxx = Sxx - (cp_index(j+1) - cp_index(j))*pow(xbar, 2.);
        Sxy = Sxy - (cp_index(j+1) - cp_index(j))*double(xbar)*ybar;
        
        
        for(int i=cp_index(j); i<cp_index(j+1); i++){
            out(i) = ybar + Sxy/Sxx * (double(i)-xbar);
        }
    }
    
    Sxx=0;
    Sxy=0;
    xbar=0;
    ybar=0;
    
    for(int i=cp_index(J-1); i<n; i++){
        xbar += i;
        ybar += data(i);
        Sxx += double(pow(i, 2.));
        Sxy += double(i)*data(i);
    }
    xbar = double(xbar) / (n - cp_index(J-1));
    ybar = ybar / (n - cp_index(J-1));
    Sxx = Sxx - (n - cp_index(J-1))*pow(double(xbar), 2.);
    Sxy = Sxy - (n - cp_index(J-1))*double(xbar)*ybar;
    
    
    for(int i=cp_index(J-1); i<n; i++){
        out(i) = ybar + Sxy/Sxx * (double(i)-xbar);
    }

    return out;
    
}

// [[Rcpp::export]]
double mosumCwald_rss(NumericVector data, NumericVector cp_index) {
    int n=data.size();
    double out=0;
    // NumericVector fit(n);
    
    // int J=cp_index.size();
    
    NumericVector fit = mosumCwald_fit(data, cp_index);
    
    for(int i=0; i<n; i++){
        out += pow(data(i)-fit(i), 2.);
    }
    
    return out;
    
}

// [[Rcpp::export]]
double slr_rss(NumericVector data) {
    int n=data.size();
    double out=0;
    NumericVector fit(n);
    
    double n_flt = double(n);
    
    double Sxx=0;
    double Sxy=0;
    double xbar=0;
    double ybar=0;
    
    for(int i=0; i<n; i++){
        xbar += i;
        ybar += data(i);
        Sxx += double(pow(i, 2.));
        Sxy += double(i)*data(i);
    }
    xbar = double(xbar) / n_flt;
    ybar = ybar / n_flt;
    Sxx = Sxx - n_flt*pow(double(xbar), 2.);
    Sxy = Sxy - n_flt*double(xbar)*ybar;
    
    
    for(int i=0; i<n; i++){
        fit(i) = ybar + Sxy/Sxx * (double(i)-xbar);
    }
    
    for(int i=0; i<n; i++){
        out += pow(data(i)-fit(i), 2.);
    }
    
    return out;
    
}

// [[Rcpp::export]]
NumericVector mosumCwald_nonsigma(NumericVector data, NumericVector time, int G) {
  int n=data.size();
  NumericVector out(n-2*G+1);
  double dt = time(n-1) / n;
  double G_flt = double(G);
  double G_half = (G_flt+1.)/2.;
  
  arma::vec diff_data(n-G);
  
  for(int i=0; i<n-G; i++){
    diff_data(i) = data(G+i) - data(i);
  }
  
  double denom = 0;
  double beta1plus = 0;
  double beta1minus = 0;
  double beta0diff = 0;
  
  /* Initializa */
  for(int i=0; i<G; i++){
    double i_flt = double(i+1);
    double w = i_flt - G_half;
    denom += pow(w, 2.0);
    beta1plus += w*data(G+i);
    beta1minus += w*data(i);
    beta0diff += diff_data(i);
  }
  
  denom *= dt;
  
  beta1plus /= denom;
  beta1minus /= denom;
  beta0diff = beta0diff/G_flt - (G_half*beta1plus + (G_half-1.)*beta1minus)*dt;
  
  double sum_plus = 0; double sum_minus = 0;
  // double sumsq_plus = 0; double sumsq_minus = 0;
  for(int i=0; i<G; i++){
    sum_plus += data(G+i);
    sum_minus += data(i);
    // sumsq_plus += pow(data(G+i), 2.);
    // sumsq_minus += pow(data(i), 2.);
  }
  
  // double SSTplus = sumsq_plus - pow(sum_plus, 2.)/G_flt;
  // double SSTminus = sumsq_minus - pow(sum_minus, 2.)/G_flt;
  // double sigmasq_plus = (SSTplus - pow(beta1plus, 2.)*denom*dt)/(G_flt-2.);
  // double sigmasq_minus = (SSTminus - pow(beta1minus, 2.)*denom*dt)/(G_flt-2.);
  // 
  // double sigma_sq = (sigmasq_plus + sigmasq_minus)/2.;
  
  out(0) = sqrt(G_flt/(8.)*(pow(beta0diff, 2.) + pow(G_flt*(beta1plus - beta1minus)*dt, 2.)/3.));
  
  /* update */
  double Delta_b1plus = 0;
  double Delta_b1minus = 0;
  double Delta_b0diff = 0;
  for(int k=G-1; k<n-G-1; k++){
    Delta_b1plus = (G_flt*data(k+G+1) - sum_plus - G_half*diff_data(k+1))/denom;
    Delta_b1minus = (G_flt*data(k+1) - sum_minus - G_half*diff_data(k-G+1))/denom;
    Delta_b0diff = (diff_data(k+1) - diff_data(k-G+1))/G_flt - (G_half*Delta_b1plus + (G_half-1.)*Delta_b1minus)*dt;
    beta1plus += Delta_b1plus;
    beta1minus += Delta_b1minus;
    beta0diff += Delta_b0diff;
    sum_plus += data(k+G+1) - data(k+1);
    sum_minus += data(k+1) - data(k-G+1);
    // sumsq_plus += pow(data(k+G+1), 2.) - pow(data(k+1), 2.);
    // sumsq_minus += pow(data(k+1), 2.) - pow(data(k-G+1), 2.);
    // SSTplus = sumsq_plus - pow(sum_plus, 2.)/G_flt;
    // SSTminus = sumsq_minus - pow(sum_minus, 2.)/G_flt;
    // sigmasq_plus = (SSTplus - pow(beta1plus, 2.)*denom*dt)/(G_flt-2.);
    // sigmasq_minus = (SSTminus - pow(beta1minus, 2.)*denom*dt)/(G_flt-2.);
    // sigma_sq = (sigmasq_plus + sigmasq_minus)/2.;
    out(k-G+2) = sqrt(G_flt/8.*(pow(beta0diff, 2.) + pow(G_flt*(beta1plus - beta1minus)*dt, 2.)/3.));
  }
  
  return out;
  
}


