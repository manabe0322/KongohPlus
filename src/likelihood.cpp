#include "header.h"

// [[Rcpp::export]]
std::vector<double> calc_eph_dedup_1_cpp(double tph,
                                         std::vector<double> mp,
                                         std::vector<double> deg,
                                         int hnc,
                                         int n_dedup_1,
                                         std::vector<double> eph_base,
                                         std::vector<std::vector<int>> index_dedup_1,
                                         std::vector<double> size_dedup_1){
  std::vector<double> eph_dedup_1(n_dedup_1, 0);
  for(int i = 0; i < hnc; ++i){
    double tph_c1 = tph * mp[i];
    double deg_c1 = deg[i];
    std::vector<int> ind_dedup1_al1 = index_dedup_1.at(2 * i);
    std::vector<int> ind_dedup1_al2 = index_dedup_1.at(2 * i + 1);
    for(int j = 0; j < n_dedup_1; ++j){
      eph_dedup_1[j] += eph_base[ind_dedup1_al1[j]] * tph_c1 * exp(deg_c1 * size_dedup_1[j]);
      eph_dedup_1[j] += eph_base[ind_dedup1_al2[j]] * tph_c1 * exp(deg_c1 * size_dedup_1[j]);
    }
  }
  return(eph_dedup_1);
}
