#include "header.h"

// [[Rcpp::export]]
std::vector<double> make_detectable_peak(std::vector<double> gtc,
                                         std::vector<double> peaks,
                                         std::vector<double> bsrs,
                                         std::vector<double> fsrs,
                                         std::vector<double> dsrs,
                                         std::vector<double> m2srs){
  std::vector<double> detectable_peak = gtc;
  int n_al = gtc.size();
  for(int i = 0; i < n_al; ++i){
    double al = gtc[i];
    int index_al = find_index_double(al, peaks);
    if(index_al != -1){
      if(bsrs[index_al] != 0){
        double bs = calc_bs_one(al);
        detectable_peak.push_back(bs);
      }
      if(fsrs[index_al] != 0){
        double fs = calc_fs_one(al);
        detectable_peak.push_back(fs);
      }
      if(dsrs[index_al] != 0){
        double ds = calc_ds_one(al);
        detectable_peak.push_back(ds);
      }
      if(m2srs[index_al] != 0){
        double m2s = calc_m2s_one(al);
        detectable_peak.push_back(m2s);
      }
    }
  }

  std::sort(detectable_peak.begin(), detectable_peak.end());
  decltype(detectable_peak)::iterator result = std::unique(detectable_peak.begin(), detectable_peak.end());
  detectable_peak.erase(result, detectable_peak.end());

  return(detectable_peak);
}

// [[Rcpp::export]]
std::vector<int> search_index_dropin(std::vector<double> peaks_rm99,
                                     std::vector<double> detectable_peak){
  int n_peaks_rm99 = peaks_rm99.size();
  std::vector<int> index_dropin;
  for(int i = 0; i < n_peaks_rm99; ++i){
    int index_peak = find_index_double(peaks_rm99[i], detectable_peak);
    if(index_peak == -1){
      index_dropin.push_back(i);
    }
  }
  return(index_dropin);
}

// [[Rcpp::export]]
std::vector<bool> judge_gtc_dropin(std::vector<double> peaks_rm99,
                                   std::vector<double> heights_rm99,
                                   std::vector<double> detectable_peak,
                                   std::vector<int> index_dropin,
                                   int n_dropin_fltr,
                                   double h_dropin_fltr,
                                   std::vector<bool> bool_gtc_di){
  std::vector<bool> bool_gtc_di_new = bool_gtc_di;
  int n_dropin = index_dropin.size();
  if(n_dropin <= n_dropin_fltr){
    for(int i = 0; i < n_dropin; ++i){
      if(heights_rm99[index_dropin[i]] > h_dropin_fltr){
        bool_gtc_di_new = bool_gtc_di;
        bool_gtc_di_new[0] = false;
        break;
      }else{
        bool_gtc_di_new[index_dropin[i] + 1] = true;
      }
    }
  }else{
    bool_gtc_di_new[0] = false;
  }

  return(bool_gtc_di_new);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> remove_dropin(std::vector<double> peaks_rm99,
                                               std::vector<double> heights_rm99,
                                               std::vector<int> index_dropin){
  int n_dropin = index_dropin.size();
  std::vector<std::vector<double>> result(2);
  for(int i = 0; i < n_dropin; ++i){
    int index_rm = index_dropin[i] - i;
    peaks_rm99.erase(peaks_rm99.begin() + index_rm);
    heights_rm99.erase(heights_rm99.begin() + index_rm);
  }
  result.at(0) = peaks_rm99;
  result.at(1) = heights_rm99;
  return(result);
}

// [[Rcpp::export]]
std::vector<bool> judge_gtc_sr(std::vector<double> gtc,
                               std::vector<double> peaks,
                               std::vector<double> peaks_rm99_rmdi,
                               std::vector<double> heights_rm99_rmdi,
                               std::vector<double> bsrs,
                               std::vector<double> fsrs,
                               std::vector<double> dsrs,
                               std::vector<double> m2srs,
                               double bsr_fltr,
                               double fsr_fltr,
                               double dsr_fltr,
                               double m2sr_fltr,
                               int n_dropin_fltr,
                               double h_dropin_fltr,
                               int n_dropin,
                               std::vector<bool> bool_gtc_di){
  std::vector<bool> bool_gtc_di_new = bool_gtc_di;
  std::vector<double> cand_stutters = setdiff_double(peaks_rm99_rmdi, gtc);
  int n_cand_stutters = cand_stutters.size();

  if(n_cand_stutters != 0){
    for(int i = 0; i < n_cand_stutters; ++i){
      double cand_stutter = cand_stutters[i];
      int index_cand_stutter = find_index_double(cand_stutter, peaks_rm99_rmdi);
      double height_stutter = heights_rm99_rmdi[index_cand_stutter];

      bool over_bsr = true;
      bool over_fsr = true;
      bool over_dsr = true;
      bool over_m2sr = true;

      double parent_allele;
      int index_parent_allele;
      int index_bsr;
      int index_fsr;
      int index_dsr;
      int index_m2sr;
      int index_di;
      double height_allele;

      parent_allele = calc_parent_al_bs_one(cand_stutter);
      index_bsr = find_index_double(parent_allele, peaks);
      if(index_bsr != -1){
        if(bsrs[index_bsr] != 0){
          index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
          if(index_parent_allele != -1){
            height_allele = heights_rm99_rmdi[index_parent_allele];
            over_bsr = height_stutter / height_allele > bsr_fltr;
          }
        }
      }

      if(over_bsr){
        parent_allele = calc_parent_al_fs_one(cand_stutter);
        index_fsr = find_index_double(parent_allele, peaks);
        if(index_fsr != -1){
          if(fsrs[index_fsr] != 0){
            index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
            if(index_parent_allele != -1){
              height_allele = heights_rm99_rmdi[index_parent_allele];
              over_fsr = height_stutter / height_allele > fsr_fltr;
            }
          }
        }

        if(over_fsr){
          parent_allele = calc_parent_al_ds_one(cand_stutter);
          index_dsr = find_index_double(parent_allele, peaks);
          if(index_dsr != -1){
            if(dsrs[index_dsr] != 0){
              index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
              if(index_parent_allele != -1){
                height_allele = heights_rm99_rmdi[index_parent_allele];
                over_dsr = height_stutter / height_allele > dsr_fltr;
              }
            }
          }

          if(over_dsr){
            parent_allele = calc_parent_al_m2s_one(cand_stutter);
            index_m2sr = find_index_double(parent_allele, peaks);
            if(index_m2sr != -1){
              if(m2srs[index_m2sr] != 0){
                index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
                if(index_parent_allele != -1){
                  height_allele = heights_rm99_rmdi[index_parent_allele];
                  over_m2sr = height_stutter / height_allele > m2sr_fltr;
                }
              }
            }

            if(over_m2sr){
              if(height_stutter <= h_dropin_fltr && n_dropin < n_dropin_fltr){
                index_di = find_index_double(cand_stutter, peaks);
                bool_gtc_di_new[index_di + 1] = true;
                n_dropin += 1;
              }else{
                bool_gtc_di_new = bool_gtc_di;
                bool_gtc_di_new[0] = false;
                break;
              }
            }
          }
        }
      }
    }
  }

  return(bool_gtc_di_new);
}

// [[Rcpp::export]]
std::vector<bool> judge_gtc(std::vector<double> gtc,
                            std::vector<double> peaks,
                            std::vector<double> heights,
                            std::vector<double> peaks_rm99,
                            std::vector<double> heights_rm99,
                            std::vector<double> bsrs,
                            std::vector<double> fsrs,
                            std::vector<double> dsrs,
                            std::vector<double> m2srs,
                            double bsr_fltr,
                            double fsr_fltr,
                            double dsr_fltr,
                            double m2sr_fltr,
                            int n_dropin_fltr,
                            double h_dropin_fltr){
  int n_peak = peaks.size();
  std::vector<bool> bool_gtc_di(n_peak + 1, false);
  bool_gtc_di[0] = true;

  std::vector<double> detectable_peak = make_detectable_peak(gtc, peaks, bsrs, fsrs, dsrs, m2srs);
  std::vector<int> index_dropin = search_index_dropin(peaks_rm99, detectable_peak);

  int n_dropin = index_dropin.size();
  if(n_dropin != 0){
    bool_gtc_di = judge_gtc_dropin(peaks_rm99, heights_rm99, detectable_peak, index_dropin, n_dropin_fltr, h_dropin_fltr, bool_gtc_di);
  }

  if(bool_gtc_di[0]){
    // Remove drop-in peaks
    std::vector<std::vector<double>> tmp = remove_dropin(peaks_rm99, heights_rm99, index_dropin);
    std::vector<double> peaks_rm99_rmdi = tmp.at(0);
    std::vector<double> heights_rm99_rmdi = tmp.at(1);

    // Judge a genotype combination considering both stutters and drop-in
    bool_gtc_di = judge_gtc_sr(gtc, peaks, peaks_rm99_rmdi, heights_rm99_rmdi,
                               bsrs, fsrs, dsrs, m2srs,
                               bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr,
                               n_dropin_fltr, h_dropin_fltr, n_dropin, bool_gtc_di);
  }

  return(bool_gtc_di);
}

// [[Rcpp::export]]
Rcpp::List judge_gtc_mk(std::vector<std::vector<double>> cand_gts,
                        std::vector<std::vector<int>> index_gtc_cpp,
                        std::vector<double> peaks,
                        std::vector<std::vector<double>> heights_all_rep,
                        std::vector<double> bsrs,
                        std::vector<double> fsrs,
                        std::vector<double> dsrs,
                        std::vector<double> m2srs,
                        double bsr_fltr,
                        double fsr_fltr,
                        double dsr_fltr,
                        double m2sr_fltr,
                        int n_dropin_fltr,
                        double h_dropin_fltr){
  int n_gt_comb = index_gtc_cpp.size();
  int hnc = index_gtc_cpp.at(0).size();
  int n_peak = peaks.size();
  int n_rep = heights_all_rep.size();

  // Remove allele Q (i.e., 99) from the object 'peaks'
  std::vector<double> peaks_rm99 = peaks;
  peaks_rm99.erase(peaks_rm99.begin() + n_peak - 1);

  std::vector<std::vector<int>> index_gtc_adopt_cpp;
  std::vector<std::vector<std::vector<bool>>> bool_gtc_di_mk_rep(n_rep);

  int count = 0;
  for(int i = 0; i < n_gt_comb; ++i){
    std::vector<int> ind_gtc = index_gtc_cpp.at(i);
    std::vector<double> gtc = extract_gtc(cand_gts, ind_gtc, hnc);

    bool adopt_gtc = true;
    std::vector<std::vector<bool>> bool_gtc_di_rep(n_rep, std::vector<bool>(n_peak + 1, false));
    for(int j = 0; j < n_rep; ++j){
      std::vector<double> heights = heights_all_rep.at(j);

      // Remove the peak height of the allele Q (i.e., 99) from the object 'heights'
      std::vector<double> heights_rm99 = heights;
      heights_rm99.erase(heights_rm99.begin() + n_peak - 1);

      bool_gtc_di_rep[j] = judge_gtc(gtc, peaks, heights, peaks_rm99, heights_rm99,
                                     bsrs, fsrs, dsrs, m2srs,
                                     bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr);
      adopt_gtc = bool_gtc_di_rep[j][0];
      if(!adopt_gtc){
        break;
      }
    }

    if(adopt_gtc){
      index_gtc_adopt_cpp.push_back(std::vector<int>(hnc));
      for(int j = 0; j < hnc; ++j){
        index_gtc_adopt_cpp[count][j] = ind_gtc[j];
      }

      for(int j = 0; j < n_rep; ++j){
        bool_gtc_di_mk_rep[j].push_back(std::vector<bool>(n_peak + 1));
        for(int k = 0; k < n_peak + 1; ++k){
          bool_gtc_di_mk_rep[j][count][k] = bool_gtc_di_rep[j][k];
        }
      }

      count += 1;
    }
  }

  List result = List::create(index_gtc_adopt_cpp, bool_gtc_di_mk_rep);
  return(result);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> extract_gtc_adopt(std::vector<std::vector<double>> cand_gts,
                                                   std::vector<std::vector<int>> index_gtc_adopt_cpp){
  int n_gtc_adopt = index_gtc_adopt_cpp.size();
  int hnc = index_gtc_adopt_cpp.at(0).size();
  std::vector<std::vector<double>> gtc_adopt(n_gtc_adopt, std::vector<double>(2 * hnc));

  for(int i = 0; i < n_gtc_adopt; ++i){
    std::vector<int> ind_gtc = index_gtc_adopt_cpp.at(i);
    for(int j = 0; j < hnc; ++j){
      std::vector<double> gt = cand_gts.at(ind_gtc[j]);
      gtc_adopt[i][2 * j] = gt[0];
      gtc_adopt[i][2 * j + 1] = gt[1];
    }
  }

  return(gtc_adopt);
}
