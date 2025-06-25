#include "header.h"

// [[Rcpp::export]]
std::vector<std::vector<int>> connect_rph_gtc(std::vector<double> gtc,
                                              std::vector<double> bs_gtc,
                                              std::vector<double> fs_gtc,
                                              std::vector<double> ds_gtc,
                                              std::vector<double> m2s_gtc,
                                              std::vector<double> peaks,
                                              std::vector<double> bsrs,
                                              std::vector<double> fsrs,
                                              std::vector<double> dsrs,
                                              std::vector<double> m2srs,
                                              std::vector<bool> bool_gtc_di){
  int n_allele = gtc.size();
  int hnc = n_allele / 2;
  int n_peak = peaks.size();
  std::vector<std::vector<int>> index_rph_gtc_cpp(n_peak, std::vector<int>(n_allele, 5 * n_peak));

  for(int i = 0; i < n_allele; ++i){
    double al = gtc[i];
    double bs = bs_gtc[i];
    double fs = fs_gtc[i];
    double ds = ds_gtc[i];
    double m2s = m2s_gtc[i];

    int ind_al = find_index_double(al, peaks);
    index_rph_gtc_cpp[ind_al][i] = ind_al;

    int ind_bs = find_index_double(bs, peaks);
    if(ind_bs != -1){
      if(bsrs[ind_al] != 0){
        index_rph_gtc_cpp[ind_bs][i] = n_peak + ind_bs;
      }
    }

    int ind_fs = find_index_double(fs, peaks);
    if(ind_fs != -1){
      if(fsrs[ind_al] != 0){
        index_rph_gtc_cpp[ind_fs][i] = 2 * n_peak + ind_fs;
      }
    }

    int ind_ds = find_index_double(ds, peaks);
    if(ind_ds != -1){
      if(dsrs[ind_al] != 0){
        index_rph_gtc_cpp[ind_ds][i] = 3 * n_peak + ind_ds;
      }
    }

    int ind_m2s = find_index_double(m2s, peaks);
    if(ind_m2s != -1){
      if(m2srs[ind_al] != 0){
        index_rph_gtc_cpp[ind_m2s][i] = 4 * n_peak + ind_m2s;
      }
    }
  }

  //Index for drop-in
  for(int i = 0; i < n_peak; ++i){
    if(bool_gtc_di[i + 1]){
      index_rph_gtc_cpp.at(i) = std::vector<int>(n_allele, 5 * n_peak + 1);
    }
  }

  //Ascending order of index for one contributor
    //When gtc1 is (12, 15), indices for peak 15 are (5 * n_peak, ind_al).
    //When gtc2 is (15, 18), indices for peak 15 are (ind_al, 5 * n_peak).
    //(5 * n_peak, ind_al) and (ind_al, 5 * n_peak) are substantially the same.
    //In that case, put (5 * n_peak, ind_al) in ascending order
  for(int i = 0; i < n_peak; ++i){
    for(int j = 0; j < hnc; ++j){
      int index_1 = index_rph_gtc_cpp[i][2 * j];
      int index_2 = index_rph_gtc_cpp[i][2 * j + 1];
      if(index_1 > index_2){
        index_rph_gtc_cpp[i][2 * j] = index_2;
        index_rph_gtc_cpp[i][2 * j + 1] = index_1;
      }
    }
  }

  return(index_rph_gtc_cpp);
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<int>>> connect_rph_gtc_mk(std::vector<std::vector<double>> cand_gts,
                                                              std::vector<std::vector<int>> index_gt_comb_adopt,
                                                              std::vector<double> peaks,
                                                              std::vector<double> bsrs,
                                                              std::vector<double> fsrs,
                                                              std::vector<double> dsrs,
                                                              std::vector<double> m2srs,
                                                              std::vector<std::vector<bool>> bool_gtc_di_mk){
  int n_gt_comb = index_gt_comb_adopt.size();
  int hnc = index_gt_comb_adopt.at(0).size();
  int n_peak = peaks.size();

  std::vector<std::vector<std::vector<int>>> index_rph_gtc_mk(n_gt_comb, std::vector<std::vector<int>>(n_peak, std::vector<int>(2 * hnc)));

  #pragma omp parallel for
  for(int i = 0; i < n_gt_comb; ++i){
    std::vector<int> ind_gtc = index_gt_comb_adopt[i];
    std::vector<double> gtc = extract_gt_comb(cand_gts, ind_gtc, hnc);
    std::vector<bool> bool_gtc_di = bool_gtc_di_mk[i];

    std::vector<double> bs_gtc = calc_bs(gtc);
    std::vector<double> fs_gtc = calc_fs(gtc);
    std::vector<double> ds_gtc = calc_ds(gtc);
    std::vector<double> m2s_gtc = calc_m2s(gtc);

    auto result = connect_rph_gtc(gtc, bs_gtc, fs_gtc, ds_gtc, m2s_gtc,
                                  peaks, bsrs, fsrs, dsrs, m2srs, bool_gtc_di);
    index_rph_gtc_mk.at(i) = std::move(result);
  }

  return(index_rph_gtc_mk);
}
