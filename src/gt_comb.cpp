#include "header.h"

bool judge_gtc_sr(std::vector<double> gtc,
                  std::vector<double> peaks_rm99,
                  std::vector<double> heights_rm99,
                  std::vector<bool> sr_consider,
                  double bsr_fltr,
                  double fsr_fltr,
                  double dsr_fltr,
                  double m2sr_fltr){
  bool adopt_gtc = true;

  bool bsr_consider = sr_consider[0];
  bool fsr_consider = sr_consider[1];
  bool dsr_consider = sr_consider[2];
  bool m2sr_consider = sr_consider[3];

  if(bsr_consider || fsr_consider || dsr_consider || m2sr_consider){
    std::vector<double> cand_stutters = setdiff_double(peaks_rm99, gtc);
    int n_cand_stutters = cand_stutters.size();

    if(n_cand_stutters != 0){
      for(int i = 0; i < n_cand_stutters; ++i){
        double cand_stutter = cand_stutters[i];
        int index_cand_stutter = find_index_double(cand_stutter, peaks_rm99);
        double height_stutter = heights_rm99[index_cand_stutter];

        bool over_bsr = true;
        bool over_fsr = true;
        bool over_dsr = true;
        bool over_m2sr = true;

        double parent_allele;
        int index_parent_allele;
        double height_allele;

        if(bsr_consider){
          parent_allele = cand_stutter + 1;
          index_parent_allele = find_index_double(parent_allele, peaks_rm99);
          if(index_parent_allele != -1){
            height_allele = heights_rm99[index_parent_allele];
            over_bsr = height_stutter / height_allele > bsr_fltr;
          }
        }

        if(over_bsr){
          if(fsr_consider){
            parent_allele = cand_stutter - 1;
            index_parent_allele = find_index_double(parent_allele, peaks_rm99);
            if(index_parent_allele != -1){
              height_allele = heights_rm99[index_parent_allele];
              over_fsr = height_stutter / height_allele > fsr_fltr;
            }
          }

          if(over_fsr){
            if(dsr_consider){
              parent_allele = cand_stutter + 2;
              index_parent_allele = find_index_double(parent_allele, peaks_rm99);
              if(index_parent_allele != -1){
                height_allele = heights_rm99[index_parent_allele];
                over_dsr = height_stutter / height_allele > dsr_fltr;
              }
            }

            if(over_dsr){
              if(m2sr_consider){
                double residual = fmod(cand_stutter, 1);
                if(residual < 0.15){
                  parent_allele = cand_stutter + 0.2;
                }else{
                  parent_allele = cand_stutter + 0.8;
                }
                index_parent_allele = find_index_double(parent_allele, peaks_rm99);
                if(index_parent_allele != -1){
                  height_allele = heights_rm99[index_parent_allele];
                  over_m2sr = height_stutter / height_allele > m2sr_fltr;
                }
              }
            }
          }
        }

        if(over_bsr && over_fsr && over_dsr && over_m2sr){
          adopt_gtc = false;
          break;
        }
      }
    }
  }
  return(adopt_gtc);
}

bool judge_gtc_phr(std::vector<double> gtc,
                   std::vector<double> peaks,
                   std::vector<double> heights,
                   double phr_fltr,
                   double st){
  bool adopt_gtc = true;
  int hnc = gtc.size() / 2;
  std::vector<int> index_gtc_all(2 * hnc);
  std::iota(index_gtc_all.begin(), index_gtc_all.end(), 0);

  for(int i = 0; i < hnc; ++i){
    double al_1 = gtc[2 * i];
    double al_2 = gtc[2 * i + 1];

    if(al_1 != al_2){
      std::vector<int> index_gtc_current(2);
      index_gtc_current[0] = 2 * i;
      index_gtc_current[1] = 2 * i + 1;
      std::vector<int> index_gtc_other = setdiff_int(index_gtc_all, index_gtc_current);

      int n_gtc_other = 2 * hnc - 2;
      std::vector<double> gtc_other(n_gtc_other);
      for(int j = 0; j < n_gtc_other; ++j){
        gtc_other[j] = gtc[index_gtc_other[j]];
      }

      int index_gtc_other_1 = find_index_double(al_1, gtc_other);
      int index_gtc_other_2 = find_index_double(al_2, gtc_other);
      if(index_gtc_other_1 == -1 && index_gtc_other_2 == -1){
        int index_peak_1 = find_index_double(al_1, peaks);
        int index_peak_2 = find_index_double(al_2, peaks);

        double h1 = heights[index_peak_1];
        double h2 = heights[index_peak_2];
        if(h1 >= st && h2 >= st){
          double phr = h1 / h2;
          if(phr < phr_fltr || phr > 1 / phr_fltr){
            adopt_gtc = false;
            break;
          }
        }
      }
    }
  }
  return(adopt_gtc);
}

// [[Rcpp::export]]
bool judge_gtc_mr(std::vector<double> gtc,
                  std::vector<double> peaks,
                  std::vector<double> heights,
                  double mr_fltr,
                  double st){
  bool adopt_gtc = true;
  int hnc = gtc.size() / 2;
  std::vector<int> index_gtc_all(2 * hnc);
  std::iota(index_gtc_all.begin(), index_gtc_all.end(), 0);
  std::vector<double> uniq_heights;
  std::vector<bool> pass_st;

  for(int i = 0; i < hnc; ++i){
    double al_1 = gtc[2 * i];
    double al_2 = gtc[2 * i + 1];

    std::vector<int> index_gtc_current(2);
    index_gtc_current[0] = 2 * i;
    index_gtc_current[1] = 2 * i + 1;
    std::vector<int> index_gtc_other = setdiff_int(index_gtc_all, index_gtc_current);

    int n_gtc_other = 2 * hnc - 2;
    std::vector<double> gtc_other(n_gtc_other);
    for(int j = 0; j < n_gtc_other; ++j){
      gtc_other[j] = gtc[index_gtc_other[j]];
    }

    int index_gtc_other_1 = find_index_double(al_1, gtc_other);
    int index_gtc_other_2 = find_index_double(al_2, gtc_other);

    double n_divide = 0;
    std::vector<double> uniq_h_per_al(2, 0);
    if(index_gtc_other_1 == -1){
      int index_peak_1 = find_index_double(al_1, peaks);
      uniq_h_per_al[0] = heights[index_peak_1];
      n_divide += 1;
    }
    if(index_gtc_other_2 == -1){
      if(al_1 != al_2){
        int index_peak_2 = find_index_double(al_2, peaks);
        uniq_h_per_al[1] = heights[index_peak_2];
      }
      n_divide += 1;
    }

    if(n_divide != 0){
      double uniq_h = accumulate(uniq_h_per_al.begin(), uniq_h_per_al.end(), 0.0) / n_divide;
      uniq_heights.push_back(uniq_h);
      pass_st.push_back(uniq_h >= st);
    }
  }

  int n_uniq = uniq_heights.size();
  if(n_uniq >= 2){
    for(int i = 0; i < (n_uniq - 1); ++i){
      if(pass_st[i]){
        double uniq_heights_target = uniq_heights[i];
        std::vector<double> uniq_heights_others = extract_subvec_double(uniq_heights, i + 1, n_uniq);
        int n_others = uniq_heights_others.size();
        for(int j = 0; j < n_others; ++j){
          double mr = uniq_heights_target / uniq_heights_others[j];
          if(mr < mr_fltr){
            adopt_gtc = false;
            break;
          }
        }
      }
    }
  }
  return(adopt_gtc);
}

// [[Rcpp::export]]
bool judge_gtc(std::vector<double> gtc,
               std::vector<double> peaks,
               std::vector<double> heights,
               std::vector<double> peaks_rm99,
               std::vector<double> heights_rm99,
               std::vector<bool> sr_consider,
               double bsr_fltr,
               double fsr_fltr,
               double dsr_fltr,
               double m2sr_fltr,
               double phr_fltr,
               double mr_fltr,
               double st){
  std::vector<double> cand_peak = make_cand_peak(gtc, sr_consider);
  bool adopt_gtc = is_element_all_vec_double(peaks_rm99, cand_peak);

  // Write other conditions
  if(adopt_gtc){
    adopt_gtc = judge_gtc_sr(gtc, peaks_rm99, heights_rm99, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr);
    if(adopt_gtc){
      adopt_gtc = judge_gtc_phr(gtc, peaks, heights, phr_fltr, st);
//      if(adopt_gtc){
//        adopt_gtc = judge_gtc_mr(gtc, peaks, heights, mr_fltr, st);
//      }
    }
  }

  return(adopt_gtc);
}

// [[Rcpp::export]]
std::vector<std::vector<int>> judge_gtc_mk(std::vector<std::vector<double>> cand_gts,
                                           std::vector<std::vector<int>> index_gtc_cpp,
                                           std::vector<double> peaks,
                                           std::vector<std::vector<double>> heights_all_rep,
                                           std::vector<bool> sr_consider,
                                           double bsr_fltr,
                                           double fsr_fltr,
                                           double dsr_fltr,
                                           double m2sr_fltr,
                                           double phr_fltr,
                                           double mr_fltr,
                                           double st){
  int n_gt_comb = index_gtc_cpp.size();
  int hnc = index_gtc_cpp.at(0).size();
  int n_peak = peaks.size();
  int n_rep = heights_all_rep.size();

  std::vector<double> peaks_rm99 = peaks;
  peaks_rm99.erase(peaks_rm99.begin() + n_peak - 1);

  std::vector<std::vector<double>> heights_all_rep_rm99(n_rep);
  for(int i = 0; i < n_rep; ++i){
    std::vector<double> heights_rm99 = heights_all_rep.at(i);
    heights_rm99.erase(heights_rm99.begin() + n_peak - 1);
    heights_all_rep_rm99.at(i) = heights_rm99;
  }

  std::vector<std::vector<int>> index_gtc_adopt_cpp;

  int count = 0;
  for(int i = 0; i < n_gt_comb; ++i){
    std::vector<int> ind_gtc = index_gtc_cpp.at(i);
    std::vector<double> gtc(2 * hnc);
    for(int j = 0; j < hnc; ++j){
      std::vector<double> gt = cand_gts.at(ind_gtc[j]);
      gtc[2 * j] = gt[0];
      gtc[2 * j + 1] = gt[1];
    }

    bool adopt_gtc = true;
    for(int j = 0; j < n_rep; ++j){
      std::vector<double> heights = heights_all_rep.at(j);
      std::vector<double> heights_rm99 = heights_all_rep_rm99.at(j);
      adopt_gtc = judge_gtc(gtc, peaks, heights, peaks_rm99, heights_rm99, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, phr_fltr, mr_fltr, st);
      if(!adopt_gtc){
        break;
      }
    }

    if(adopt_gtc){
      index_gtc_adopt_cpp.push_back(std::vector<int>(hnc));
      for(int k = 0; k < hnc; ++k){
        index_gtc_adopt_cpp[count][k] = ind_gtc[k];
      }
      count += 1;
    }
  }

  return(index_gtc_adopt_cpp);
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
