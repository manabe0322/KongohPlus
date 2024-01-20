#include "header.h"

std::vector<double> calc_eph_base_gt1(std::vector<double> gt,
                                      std::vector<double> peaks,
                                      std::vector<double> bs_all,
                                      std::vector<double> fs_all,
                                      std::vector<double> ds_all,
                                      std::vector<double> m2s_all,
                                      std::vector<double> bsrs,
                                      std::vector<double> fsrs,
                                      std::vector<double> dsrs,
                                      double m2sr,
                                      double ae){

  int n_peak = peaks.size();
  int n_al = gt.size();

  std::vector<double> eph_base_gt1(n_peak);

  for(int i = 0; i < n_al; ++i){
    double al = gt[i];
    int index_al = find_index_double(al, peaks);
    double bsr = bsrs[index_al];
    double fsr = fsrs[index_al];
    double dsr = dsrs[index_al];

    eph_base_gt1[index_al] += ae / (1 + bsr + fsr + dsr + m2sr);

    double bs = bs_all[index_al];
    int index_bs = find_index_double(bs, peaks);
    if(index_bs != -1){
      eph_base_gt1[index_bs] += ae * bsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double fs = fs_all[index_al];
    int index_fs = find_index_double(fs, peaks);
    if(index_fs != -1){
      eph_base_gt1[index_fs] += ae * fsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double ds = ds_all[index_al];
    int index_ds = find_index_double(ds, peaks);
    if(index_ds != -1){
      eph_base_gt1[index_ds] += ae * dsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double m2s = m2s_all[index_al];
    int index_m2s = find_index_double(m2s, peaks);
    if(index_m2s != -1){
      eph_base_gt1[index_m2s] += ae * m2sr / (1 + bsr + fsr + dsr + m2sr);
    }
  }
  return(eph_base_gt1);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> make_eph_base(std::vector<double> peaks,
                                               std::vector<std::vector<double>> cand_gts,
                                               std::vector<bool> sr_consider,
                                               std::vector<double> bsrs,
                                               std::vector<double> fsrs,
                                               std::vector<double> dsrs,
                                               double m2sr,
                                               double ae){
  int n_cand_gts = cand_gts.size();
  int n_peaks = peaks.size();

  std::vector<double> bs_all = calc_bs(peaks);
  std::vector<double> fs_all = calc_fs(peaks);
  std::vector<double> ds_all = calc_ds(peaks);
  std::vector<double> m2s_all = calc_m2s(peaks);

  std::vector<std::vector<double>> eph_base(n_cand_gts, std::vector<double>(n_peaks));

  for(int i = 0; i < n_cand_gts; ++i){
    std::vector<double> gt = cand_gts.at(i);
    eph_base.at(i) = calc_eph_base_gt1(gt, peaks, bs_all, fs_all, ds_all, m2s_all, bsrs, fsrs, dsrs, m2sr, ae);
  }

  return(eph_base);
}
