#include "header.h"

// [[Rcpp::export]]
std::vector<int> extract_subvec_int(std::vector<int> vec, int index_start, int index_end){
  std::vector<int> subvec(vec.begin() + index_start, vec.begin() + index_end + 1);
  return(subvec);
}

// [[Rcpp::export]]
std::vector<double> extract_subvec_double(std::vector<double> vec, int index_start, int index_end){
  std::vector<double> subvec(vec.begin() + index_start, vec.begin() + index_end + 1);
  return(subvec);
}

// [[Rcpp::export]]
int find_index_double(double element, std::vector<double> vec){
  int index = -1;
  auto tmp = std::find(vec.begin(), vec.end(), element);
  if(tmp != vec.end()){
    index = tmp - vec.begin();
  }
  return(index);
}

// [[Rcpp::export]]
bool is_element_all_vec_double(std::vector<double> vec_a, std::vector<double> vec_b){
  bool ans = true;
  int n_vec_a = vec_a.size();
  for(int i = 0; i < n_vec_a; ++i){
    double element = vec_a[i];
    auto index = std::find(vec_b.begin(), vec_b.end(), element);
    ans = index != vec_b.end();
    if(!ans){
      break;
    }
  }
  return(ans);
}

// [[Rcpp::export]]
std::vector<double> setdiff_double(std::vector<double> vec_num1, std::vector<double> vec_num2){
  std::vector<double> only_vec_num1;
  int n_vec1 = vec_num1.size();
  int n_vec2 = vec_num2.size();
  for(int i = 0; i < n_vec1; ++i){
    double vec1_1 = vec_num1[i];
    bool only = true;
    for(int j = 0; j < n_vec2; ++j){
      double vec2_1 = vec_num2[j];
      if(vec1_1 == vec2_1){
        only = false;
        break;
      }
    }
    if(only){
      only_vec_num1.push_back(vec1_1);
    }
  }
  return(only_vec_num1);
}

// [[Rcpp::export]]
std::vector<int> setdiff_int(std::vector<int> vec_int1, std::vector<int> vec_int2){
  std::vector<int> only_vec_int1;
  int n_vec1 = vec_int1.size();
  int n_vec2 = vec_int2.size();
  for(int i = 0; i < n_vec1; ++i){
    int vec1_1 = vec_int1[i];
    bool only = true;
    for(int j = 0; j < n_vec2; ++j){
      int vec2_1 = vec_int2[j];
      if(vec1_1 == vec2_1){
        only = false;
        break;
      }
    }
    if(only){
      only_vec_int1.push_back(vec1_1);
    }
  }
  return(only_vec_int1);
}

// [[Rcpp::export]]
std::vector<double> calc_bs(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> bs = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    bs[i] = alleles[i] - 1;
  }
  return(bs);
}

// [[Rcpp::export]]
std::vector<double> calc_fs(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> fs = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    fs[i] = alleles[i] + 1;
  }
  return(fs);
}

// [[Rcpp::export]]
std::vector<double> calc_ds(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> ds = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    ds[i] = alleles[i] - 2;
  }
  return(ds);
}

// [[Rcpp::export]]
std::vector<double> calc_m2s(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> m2s = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    double al = alleles[i];
    double residual = fmod(al, 1);
    if(residual < 0.15){
      m2s[i] = al - 0.8;
    }else{
      m2s[i] = al - 0.2;
    }
  }
  return(m2s);
}

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
