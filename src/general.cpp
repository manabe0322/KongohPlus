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
double calc_bs_one(double al){
  double tmp = al - 1;
  double bs_one = round(tmp * 10) / 10;
  return(bs_one);
}

// [[Rcpp::export]]
std::vector<double> calc_bs(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> bs = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    bs[i] = calc_bs_one(alleles[i]);
  }
  return(bs);
}

// [[Rcpp::export]]
double calc_fs_one(double al){
  double tmp = al + 1;
  double fs_one = round(tmp * 10) / 10;
  return(fs_one);
}

// [[Rcpp::export]]
std::vector<double> calc_fs(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> fs = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    fs[i] = calc_fs_one(alleles[i]);
  }
  return(fs);
}

// [[Rcpp::export]]
double calc_ds_one(double al){
  double tmp = al - 2;
  double ds_one = round(tmp * 10) / 10;
  return(ds_one);
}

// [[Rcpp::export]]
std::vector<double> calc_ds(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> ds = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    ds[i] = calc_ds_one(alleles[i]);
  }
  return(ds);
}

// [[Rcpp::export]]
double calc_m2s_one(double al){
  double residual = fmod(al, 1);
  double m2s_one = 0;
  if(residual < 0.15){
    double tmp = al - 0.8;
    m2s_one = round(tmp * 10) / 10;
  }else{
    double tmp = al - 0.2;
    m2s_one = round(tmp * 10) / 10;
  }
  return(m2s_one);
}

// [[Rcpp::export]]
std::vector<double> calc_m2s(std::vector<double> alleles){
  int n_al = alleles.size();
  std::vector<double> m2s = std::vector<double>(n_al);
  for(int i = 0; i < n_al; ++i){
    m2s[i] = calc_m2s_one(alleles[i]);
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
