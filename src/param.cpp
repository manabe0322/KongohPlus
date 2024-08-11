#include "header.h"

// [[Rcpp::export]]
std::vector<double> make_cand_peak(std::vector<double> alleles,
                                   std::vector<bool> sr_consider){
  std::vector<double> cand_peak = alleles;

  if(sr_consider[0]){
    std::vector<double> bs = calc_bs(alleles);
    int n_bs = bs.size();
    for(int i = 0; i < n_bs; ++i){
      cand_peak.push_back(bs[i]);
    }
  }

  if(sr_consider[1]){
    std::vector<double> fs = calc_fs(alleles);
    int n_fs = fs.size();
    for(int i = 0; i < n_fs; ++i){
      cand_peak.push_back(fs[i]);
    }
  }

  if(sr_consider[2]){
    std::vector<double> ds = calc_ds(alleles);
    int n_ds = ds.size();
    for(int i = 0; i < n_ds; ++i){
      cand_peak.push_back(ds[i]);
    }
  }

  if(sr_consider[3]){
    std::vector<double> m2s = calc_m2s(alleles);
    int n_m2s = m2s.size();
    for(int i = 0; i < n_m2s; ++i){
      cand_peak.push_back(m2s[i]);
    }
  }

  std::sort(cand_peak.begin(), cand_peak.end());
  decltype(cand_peak)::iterator result = std::unique(cand_peak.begin(), cand_peak.end());
  cand_peak.erase(result, cand_peak.end());

  return(cand_peak);
}
