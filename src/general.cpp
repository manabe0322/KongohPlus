#include "header.h"

// [[Rcpp::export]]
std::vector<double> calc_bs(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> bs = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    bs[i] = alleles[i] - 1;
  }
  return(bs);
}

// [[Rcpp::export]]
std::vector<double> calc_fs(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> fs = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    fs[i] = alleles[i] + 1;
  }
  return(fs);
}

// [[Rcpp::export]]
std::vector<double> calc_ds(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> ds = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    ds[i] = alleles[i] - 2;
  }
  return(ds);
}

// [[Rcpp::export]]
std::vector<double> calc_m2s(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> m2s = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    double al = alleles[i];
    double residual = fmod(al, 1);
    if(residual < 0.15){
      m2s[i] = al - 0.8;
    }else{
      m2s[i] = al - 0.2;
    }
  }
  return(m2s);
}
