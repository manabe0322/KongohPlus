#include "header.h"

bool judge_gtc_sr(std::vector<double> gtc,
                  std::vector<double> peaks_rm99,
                  std::vector<double> heights_rm99,
                  double bsr_fltr,
                  double fsr_fltr,
                  double dsr_fltr,
                  double m2sr_fltr,
                  std::vector<bool> sr_consider){
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

// [[Rcpp::export]]
bool judge_gtc(std::vector<double> gtc,
               std::vector<double> peaks_rm99,
               std::vector<double> heights_rm99,
               double bsr_fltr,
               double fsr_fltr,
               double dsr_fltr,
               double m2sr_fltr,
               std::vector<bool> sr_consider){
  std::vector<double> cand_peak = make_cand_peak(gtc, sr_consider);
  bool adopt_gtc = is_element_all_vec_double(peaks_rm99, cand_peak);

  // Write other conditions
  if(adopt_gtc){
    adopt_gtc = judge_gtc_sr(gtc, peaks_rm99, heights_rm99, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, sr_consider);
  }

  return(adopt_gtc);
}

// [[Rcpp::export]]
std::vector<std::vector<int>> judge_gt_comb_mk(std::vector<std::vector<double>> cand_gts,
                                               std::vector<std::vector<int>> index_gt_comb,
                                               std::vector<double> peaks,
                                               std::vector<double> heights,
                                               double bsr_fltr,
                                               double fsr_fltr,
                                               double dsr_fltr,
                                               double m2sr_fltr,
                                               std::vector<bool> sr_consider){
  int n_gt_comb = index_gt_comb.size();
  int hnc = index_gt_comb.at(0).size();
  int n_peak = peaks.size();

  std::vector<double> peaks_rm99 = peaks;
  peaks_rm99.erase(peaks_rm99.begin() + n_peak - 1);

  std::vector<double> heights_rm99 = heights;
  heights_rm99.erase(heights_rm99.begin() + n_peak - 1);

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

    bool adopt_gtc = judge_gtc(gtc, peaks_rm99, heights_rm99, bsr_fltr, fsr_fltr, dsr_fltr, m2sr_fltr, sr_consider);
    if(adopt_gtc){
      index_gt_comb_adopt.push_back(std::vector<int>(hnc));
      for(int k = 0; k < hnc; ++k){
        index_gt_comb_adopt[count][k] = ind_gtc[k];
      }
      count += 1;
    }
  }

  return(index_gt_comb_adopt);
}
