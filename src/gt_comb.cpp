#include "header.h"

/*
// [[Rcpp::export]]
bool judge_gtc(std::vector<double> peaks,
               std::vector<double> gtc,
               std::vector<bool> sr_consider){
  std::vector<double> cand_peak = make_cand_peak(gtc, sr_consider);
  bool adopt_gtc = is_element_all_vec_double(cand_peak, peaks);

  // Write other conditions

  return(adopt_gtc);
}

std::vector<double> calc_eph(std::vector<double> gt,
                             std::vector<double> peaks,
                             std::vector<double> bs_all,
                             std::vector<double> fs_all,
                             std::vector<double> ds_all,
                             std::vector<double> m2s_all,
                             std::vector<double> bsrs,
                             std::vector<double> fsrs,
                             std::vector<double> dsrs,
                             double m2sr){

  int n_peak = peaks.size();
  int n_al = gt.size();

  std::vector<double> eph(n_peak);

  for(int i = 0; i < n_al; ++i){
    double al = gt[i];
    int index_al = find_index_double(al, peaks);
    double bsr = bsrs[index_al];
    double fsr = fsrs[index_al];
    double dsr = dsrs[index_al];

    eph[index_al] += 1 / (1 + bsr + fsr + dsr + m2sr);

    double bs = bs_all[index_al];
    int index_bs = find_index_double(bs, peaks);
    if(index_bs != -1){
      eph[index_bs] += bsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double fs = fs_all[index_al];
    int index_fs = find_index_double(fs, peaks);
    if(index_fs != -1){
      eph[index_fs] += fsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double ds = ds_all[index_al];
    int index_ds = find_index_double(ds, peaks);
    if(index_ds != -1){
      eph[index_ds] += dsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double m2s = m2s_all[index_al];
    int index_m2s = find_index_double(m2s, peaks);
    if(index_m2s != -1){
      eph[index_m2s] += m2sr / (1 + bsr + fsr + dsr + m2sr);
    }
  }
  return(eph);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> tmp(std::vector<double> peaks,
                                     std::vector<std::vector<double>> cand_gts,
                                     std::vector<std::vector<int>> index_gt_comb,
                                     std::vector<bool> sr_consider){
  int n_gt_comb = index_gt_comb.size();
  int hnc = index_gt_comb.at(0).size();

  // std::vector<int> adopt_index;
  // adopt_index.reserve(n_gt_comb);
  for(int i = 0; i < n_gt_comb; ++i){
    std::vector<int> ind_gtc = index_gt_comb.at(i);
    std::vector<double> gtc(2 * hnc);
    for(int j = 0; j < hnc; ++j){
      std::vector<double> gt = cand_gts.at(ind_gtc[j]);
      gtc[2 * j] = gt[0];
      gtc[2 * j + 1] = gt[1];
    }

    bool adopt_gtc = judge_gtc(peaks, gtc, sr_consider);
    if(adopt_gtc){
      // adopt_index.push_back(i);
    }
  }
}
*/
