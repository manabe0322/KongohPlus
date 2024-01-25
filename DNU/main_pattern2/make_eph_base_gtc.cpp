// [[Rcpp::export]]
std::vector<std::vector<double>> make_eph_base_gtc(std::vector<std::vector<double>> eph_base,
                                                   std::vector<std::vector<int>> index_gt_comb_adopt,
                                                   int n_peak){
  int n_gtc_adopt = index_gt_comb_adopt.size();
  int hnc = index_gt_comb_adopt.at(0).size();

  std::vector<std::vector<double>> eph_base_gtc(hnc, std::vector<double>(n_gtc_adopt * n_peak));
  int count = 0;
  for(int i = 0; i < n_gtc_adopt; ++i){
    std::vector<int> ind_gtc = index_gt_comb_adopt.at(i);
    for(int j = 0; j < hnc; ++j){
      std::vector<double> eph_base_gt = eph_base.at(ind_gtc[j]);
      for(int k = 0; k < n_peak; ++k){
        eph_base_gtc.at(j).at(count + k) = eph_base_gt[k];
      }
    }
    count += n_peak;
  }

  return(eph_base_gtc);
}
