#include "header.h"

// [[Rcpp::export]]
bool judge_gtc(std::vector<double> peaks,
               std::vector<double> gtc,
               std::vector<bool> sr_consider){
  std::vector<double> cand_peak = make_cand_peak(gtc, sr_consider);
  bool adopt_gtc = is_element_all_vec_double(cand_peak, peaks);

  // Write other conditions

  return(adopt_gtc);
}

std::vector<double> calc_eph_base_c1(std::vector<double> gt,
                                     std::vector<double> peaks,
                                     std::vector<double> bs_all,
                                     std::vector<double> fs_all,
                                     std::vector<double> ds_all,
                                     std::vector<double> m2s_all,
                                     std::vector<double> bsrs,
                                     std::vector<double> fsrs,
                                     std::vector<double> dsrs,
                                     double m2sr,
                                     double ave_het_h){

  int n_peak = peaks.size();
  int n_al = gt.size();

  std::vector<double> eph_base_c1(n_peak);

  for(int i = 0; i < n_al; ++i){
    double al = gt[i];
    int index_al = find_index_double(al, peaks);
    double bsr = bsrs[index_al];
    double fsr = fsrs[index_al];
    double dsr = dsrs[index_al];

    eph_base_c1[index_al] += ave_het_h / (1 + bsr + fsr + dsr + m2sr);

    double bs = bs_all[index_al];
    int index_bs = find_index_double(bs, peaks);
    if(index_bs != -1){
      eph_base_c1[index_bs] += ave_het_h * bsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double fs = fs_all[index_al];
    int index_fs = find_index_double(fs, peaks);
    if(index_fs != -1){
      eph_base_c1[index_fs] += ave_het_h * fsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double ds = ds_all[index_al];
    int index_ds = find_index_double(ds, peaks);
    if(index_ds != -1){
      eph_base_c1[index_ds] += ave_het_h * dsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double m2s = m2s_all[index_al];
    int index_m2s = find_index_double(m2s, peaks);
    if(index_m2s != -1){
      eph_base_c1[index_m2s] += ave_het_h * m2sr / (1 + bsr + fsr + dsr + m2sr);
    }
  }
  return(eph_base_c1);
}

std::vector<std::vector<double>> calc_eph_base(std::vector<double> gtc,
                                               std::vector<double> peaks,
                                               std::vector<double> bs_all,
                                               std::vector<double> fs_all,
                                               std::vector<double> ds_all,
                                               std::vector<double> m2s_all,
                                               std::vector<double> bsrs,
                                               std::vector<double> fsrs,
                                               std::vector<double> dsrs,
                                               double m2sr,
                                               double ave_het_h){
  int hnc = gtc.size() / 2;
  std::vector<std::vector<double>> eph_base;
  for(int i = 0; i < hnc; ++i){
    std::vector<double> gt(2);
    gt[0] = gtc[2 * i];
    gt[1] = gtc[2 * i + 1];

    std::vector<double> eph_base_cl = calc_eph_base_c1(gt, peaks, bs_all, fs_all, ds_all, m2s_all, bsrs, fsrs, dsrs, m2sr, ave_het_h);
    eph_base.at(i) = eph_base_cl;
  }
  return(eph_base);
}

std::vector<double> calc_eph(std::vector<std::vector<double>> eph_base,
                             std::vector<double> mp){
  int hnc = eph_base.size();
  int n_peak = eph_base.at(0).size();
  std::vector<double> eph(n_peak);

  for(int i = 0; i < hnc; ++i){
    std::vector<double> eph_base_c1 = eph_base.at(i);
    double mp_c1 = mp[i];

    for(int j = 0; j < n_peak; ++j){
      eph[j] += eph_base_c1[j] * mp_c1;
    }
  }

  return(eph);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> tmp(std::vector<double> peaks,
                                     std::vector<std::vector<double>> cand_gts,
                                     std::vector<std::vector<int>> index_gt_comb,
                                     std::vector<bool> sr_consider,
                                     std::vector<double> bsrs,
                                     std::vector<double> fsrs,
                                     std::vector<double> dsrs,
                                     double m2sr,
                                     std::vector<std::vector<double>> mp_cand,
                                     std::vector<double> oph){
  int n_gt_comb = index_gt_comb.size();
  int hnc = index_gt_comb.at(0).size();
  int n_mp_cand = mp_cand.size();

  std::vector<double> bs_all = calc_bs(peaks);
  std::vector<double> fs_all = calc_fs(peaks);
  std::vector<double> ds_all = calc_ds(peaks);
  std::vector<double> m2s_all = calc_m2s(peaks);

  double ave_het_h = std::accumulate(oph.begin(), oph.end(), 0.0) / 2;

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
      std::vector<std::vector<double>> eph_base = calc_eph_base(gtc, peaks, bs_all, fs_all, ds_all, m2s_all, bsrs, fsrs, dsrs, m2sr, ave_het_h);

      for(int k = 0; k < n_mp_cand; ++k){
        std::vector<double> mp = mp_cand.at(k);
        std::vector<double> eph = calc_eph(eph_base, mp);
      }

      // adopt_index.push_back(i);
    }
  }
}
