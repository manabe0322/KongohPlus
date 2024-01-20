#include "header.h"

// [[Rcpp::export]]
bool judge_gtc(std::vector<double> peaks_rm99,
               std::vector<double> gtc,
               std::vector<bool> sr_consider){
  std::vector<double> cand_peak = make_cand_peak(gtc, sr_consider);
  bool adopt_gtc = is_element_all_vec_double(peaks_rm99, cand_peak);

  // Write other conditions

  return(adopt_gtc);
}

// [[Rcpp::export]]
std::vector<std::vector<int>> judge_gt_comb_mk(std::vector<double> peaks,
                                               std::vector<std::vector<double>> cand_gts,
                                               std::vector<std::vector<int>> index_gt_comb,
                                               std::vector<bool> sr_consider){
  int n_gt_comb = index_gt_comb.size();
  int hnc = index_gt_comb.at(0).size();
  int n_peak = peaks.size();

  std::vector<double> peaks_rm99 = peaks;
  peaks_rm99.erase(peaks_rm99.begin() + n_peak - 1);

  std::vector<std::vector<int>> index_gt_comb_adopt;

  int count = 0;
  for(int i = 0; i < n_gt_comb; ++i){
    std::vector<int> ind_gtc = index_gt_comb.at(i);
    std::vector<double> gtc(2 * hnc);
    for(int j = 0; j < hnc; ++j){
      std::vector<double> gt = cand_gts.at(ind_gtc[j]);
      gtc[2 * j] = gt[0];
      gtc[2 * j + 1] = gt[1];
    }

    bool adopt_gtc = judge_gtc(peaks_rm99, gtc, sr_consider);
    if(adopt_gtc){
      index_gt_comb_adopt.push_back(std::vector<int>(hnc));
      for(int k = 0; k < hnc; ++k){
        index_gt_comb_adopt[count][k] = ind_gtc[k];
      }
    }
  }

  return(index_gt_comb_adopt);
}
