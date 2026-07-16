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
double calc_parent_al_bs_one(double bs){
  double tmp = bs + 1;
  double par_al = round(tmp * 10) / 10;
  return(par_al);
}

// [[Rcpp::export]]
std::vector<double> calc_parent_al_bs(std::vector<double> bss){
  int n_bs = bss.size();
  std::vector<double> alleles = std::vector<double>(n_bs);
  for(int i = 0; i < n_bs; ++i){
    alleles[i] = calc_parent_al_bs_one(bss[i]);
  }
  return(alleles);
}

// [[Rcpp::export]]
double calc_parent_al_fs_one(double fs){
  double tmp = fs - 1;
  double par_al = round(tmp * 10) / 10;
  return(par_al);
}

// [[Rcpp::export]]
std::vector<double> calc_parent_al_fs(std::vector<double> fss){
  int n_fs = fss.size();
  std::vector<double> alleles = std::vector<double>(n_fs);
  for(int i = 0; i < n_fs; ++i){
    alleles[i] = calc_parent_al_fs_one(fss[i]);
  }
  return(alleles);
}

// [[Rcpp::export]]
double calc_parent_al_ds_one(double ds){
  double tmp = ds + 2;
  double par_al = round(tmp * 10) / 10;
  return(par_al);
}

// [[Rcpp::export]]
std::vector<double> calc_parent_al_ds(std::vector<double> dss){
  int n_ds = dss.size();
  std::vector<double> alleles = std::vector<double>(n_ds);
  for(int i = 0; i < n_ds; ++i){
    alleles[i] = calc_parent_al_ds_one(dss[i]);
  }
  return(alleles);
}

// [[Rcpp::export]]
double calc_parent_al_m2s_one(double m2s){
  double residual = fmod(m2s, 1);
  double par_al = 0;
  if(residual < 0.15){
    double tmp = m2s + 0.2;
    par_al = round(tmp * 10) / 10;
  }else{
    double tmp = m2s + 0.8;
    par_al = round(tmp * 10) / 10;
  }
  return(par_al);
}

// [[Rcpp::export]]
std::vector<double> calc_parent_al_m2s(std::vector<double> m2ss){
  int n_m2s = m2ss.size();
  std::vector<double> alleles = std::vector<double>(n_m2s);
  for(int i = 0; i < n_m2s; ++i){
    alleles[i] = calc_parent_al_m2s_one(m2ss[i]);
  }
  return(alleles);
}

// [[Rcpp::export]]
std::vector<double> extract_gt_comb(std::vector<std::vector<double>> cand_gt,
                                    std::vector<int> ind_gtc,
                                    int hnc){
  std::vector<double> gtc(2 * hnc);
  for(int i = 0; i < hnc; ++i){
    std::vector<double> gt = cand_gt[ind_gtc[i]];
    gtc[2 * i] = gt[0];
    gtc[2 * i + 1] = gt[1];
  }
  return(gtc);
}

// [[Rcpp::export]]
double log_sum_exp(const std::vector<double>& log_values){
  double max_val = *std::max_element(log_values.begin(), log_values.end());
  double sum_exp = 0.0;
  for(double v : log_values){
    sum_exp += std::exp(v - max_val);
  }
  return max_val + std::log(sum_exp);
}
