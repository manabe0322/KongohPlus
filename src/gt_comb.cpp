#include "header.h"

/*// [[Rcpp::export]]
bool judge_gt_comb(){

}*/

// [[Rcpp::export]]
std::vector<double> calc_bs_cpp(std::vector<double> al){
  int n_al = al.size();
  std::vector<double> bs = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    bs[i] = al[i] - 1;
  }
  return(bs);
}

/*
std::vector<double> make_cand_peak_cpp(std::vector<double> alleles,
                                       std::vector<bool> sr_consider){
  std::vector<double> cand_peak = alleles;

  if(sr_consider[0]){

  }

  if(sr_consider[1]){

  }

  if(sr_consider[2]){

  }

  if(sr_consider[3]){

  }
}*/

/*
std::vector<std::vector<double>> tmp(std::vector<std::vector<double>> cand_gts,
                                     std::vector<std::vector<double>> index_gt){

}*/
