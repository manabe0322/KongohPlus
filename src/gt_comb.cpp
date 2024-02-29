#include "header.h"

// [[Rcpp::export]]
std::vector<int> search_index_dropin(std::vector<double> peaks_rm99,
                                     std::vector<double> cand_peak){
  int n_peaks_rm99 = peaks_rm99.size();
  std::vector<int> index_dropin;
  for(int i = 0; i < n_peaks_rm99; ++i){
    int index_peak = find_index_double(peaks_rm99[i], cand_peak);
    if(index_peak == -1){
      index_dropin.push_back(i);
    }
  }
  return(index_dropin);
}

// [[Rcpp::export]]
bool judge_gtc_dropin(std::vector<double> peaks_rm99,
                      std::vector<double> heights_rm99,
                      std::vector<double> cand_peak,
                      std::vector<int> index_dropin,
                      int n_dropin_fltr,
                      double h_dropin_fltr){
  bool adopt_gtc = true;
  int n_dropin = index_dropin.size();
  if(n_dropin <= n_dropin_fltr){
    for(int i = 0; i < n_dropin; ++i){
      if(heights_rm99[index_dropin[i]] > h_dropin_fltr){
        adopt_gtc = false;
        break;
      }
    }
  }else{
    adopt_gtc = false;
  }
  return(adopt_gtc);
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
bool judge_gtc_sr(std::vector<double> gtc,
                  std::vector<double> peaks_rm99_rmdi,
                  std::vector<double> heights_rm99_rmdi,
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
        double height_allele;

        if(bsr_consider){
          parent_allele = cand_stutter + 1;
          index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
          if(index_parent_allele != -1){
            height_allele = heights_rm99_rmdi[index_parent_allele];
            over_bsr = height_stutter / height_allele > bsr_fltr;
          }
        }

        if(over_bsr){
          if(fsr_consider){
            parent_allele = cand_stutter - 1;
            index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
            if(index_parent_allele != -1){
              height_allele = heights_rm99_rmdi[index_parent_allele];
              over_fsr = height_stutter / height_allele > fsr_fltr;
            }
          }

          if(over_fsr){
            if(dsr_consider){
              parent_allele = cand_stutter + 2;
              index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
              if(index_parent_allele != -1){
                height_allele = heights_rm99_rmdi[index_parent_allele];
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
                index_parent_allele = find_index_double(parent_allele, peaks_rm99_rmdi);
                if(index_parent_allele != -1){
                  height_allele = heights_rm99_rmdi[index_parent_allele];
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
               int n_dropin_fltr,
               double h_dropin_fltr){
  std::vector<double> cand_peak = make_cand_peak(gtc, sr_consider);
  std::vector<int> index_dropin = search_index_dropin(peaks_rm99, cand_peak);
  bool adopt_gtc = judge_gtc_dropin(peaks_rm99, heights_rm99, cand_peak, index_dropin, n_dropin_fltr, h_dropin_fltr);

  if(adopt_gtc){
    std::vector<std::vector<double>> tmp = remove_dropin(peaks_rm99, heights_rm99, index_dropin);
    std::vector<double> peaks_rm99_rmdi = tmp.at(0);
    std::vector<double> heights_rm99_rmdi = tmp.at(1);
    adopt_gtc = judge_gtc_sr(gtc, peaks_rm99_rmdi, heights_rm99_rmdi, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr);
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
                                           int n_dropin_fltr,
                                           double h_dropin_fltr){
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
      adopt_gtc = judge_gtc(gtc, peaks, heights, peaks_rm99, heights_rm99, sr_consider, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, n_dropin_fltr, h_dropin_fltr);
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
