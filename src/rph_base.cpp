#include "header.h"

// [[Rcpp::export]]
std::vector<double> make_rph_base(std::vector<double> peaks,
                                  std::vector<double> bsrs,
                                  std::vector<double> fsrs,
                                  std::vector<double> dsrs,
                                  double m2sr){
  int n_peaks = peaks.size();
  std::vector<double> bs_all = calc_bs(peaks);
  std::vector<double> fs_all = calc_fs(peaks);
  std::vector<double> ds_all = calc_ds(peaks);
  std::vector<double> m2s_all = calc_m2s(peaks);
  std::vector<double> rph_base(5 * n_peaks + 2, 0);
  rph_base[5 * n_peaks + 1] = -1; // for drop-in

  for(int i = 0; i < n_peaks; ++i){
    double bsr = bsrs[i];
    double fsr = fsrs[i];
    double dsr = dsrs[i];

    rph_base[i] = 1 / (1 + bsr + fsr + dsr + m2sr);

    double bs = bs_all[i];
    int index_bs = find_index_double(bs, peaks);
    if(index_bs != -1){
      rph_base[n_peaks + index_bs] = bsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double fs = fs_all[i];
    int index_fs = find_index_double(fs, peaks);
    if(index_fs != -1){
      rph_base[2 * n_peaks + index_fs] = fsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double ds = ds_all[i];
    int index_ds = find_index_double(ds, peaks);
    if(index_ds != -1){
      rph_base[3 * n_peaks + index_ds] = dsr / (1 + bsr + fsr + dsr + m2sr);
    }

    double m2s = m2s_all[i];
    int index_m2s = find_index_double(m2s, peaks);
    if(index_m2s != -1){
      rph_base[4 * n_peaks + index_m2s] = m2sr / (1 + bsr + fsr + dsr + m2sr);
    }
  }
  return(rph_base);
}
