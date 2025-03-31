#include "header.h"

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_rph_dedup1_cpp(std::vector<double> mp,
                                                     std::vector<double> deg,
                                                     int hnc,
                                                     int n_rph_sum_rep1,
                                                     std::vector<double> rph_base_use,
                                                     std::vector<std::vector<int>> index_dedup1_rep1_cpp,
                                                     std::vector<double> size_rep1,
                                                     double min_size){
  std::vector<std::vector<double>> rph_dedup1(2 * hnc, std::vector<double>(n_rph_sum_rep1, 0));
  for(int i = 0; i < hnc; ++i){
    double mp_c1 = mp[i];
    double deg_c1 = deg[i];
    std::vector<int> ind_dedup1_al1 = index_dedup1_rep1_cpp.at(2 * i);
    std::vector<int> ind_dedup1_al2 = index_dedup1_rep1_cpp.at(2 * i + 1);
    for(int j = 0; j < n_rph_sum_rep1; ++j){
      double rph_mpdeg = mp_c1 * std::pow(deg_c1, (size_rep1[j] - min_size) / min_size);
      rph_dedup1[2 * i][j] = rph_base_use[ind_dedup1_al1[j]] * rph_mpdeg;
      rph_dedup1[2 * i + 1][j] = rph_base_use[ind_dedup1_al2[j]] * rph_mpdeg;
    }
  }
  return(rph_dedup1);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_gamma_param_cpp(std::vector<std::vector<double>> rph_dedup1,
                                                      std::vector<std::vector<double>> cv_square_rep1,
                                                      int n_roop,
                                                      int n_rph_sum_rep1,
                                                      double tph){
  std::vector<double> tmp1_sum(n_rph_sum_rep1, 0);
  std::vector<double> tmp2_sum(n_rph_sum_rep1, 0);
  for(int i = 0; i < n_roop; ++i){
    std::vector<double> rph_dedup1_al1 = rph_dedup1.at(i);
    std::vector<double> cv_square_rep1_al1 = cv_square_rep1.at(i);
    for(int j = 0; j < n_rph_sum_rep1; ++j){
      tmp1_sum[j] += rph_dedup1_al1[j];
      tmp2_sum[j] += rph_dedup1_al1[j] * cv_square_rep1_al1[j];
    }
  }

  std::vector<std::vector<double>> gamma_param(2, std::vector<double>(n_rph_sum_rep1, 0));
  for(int i = 0; i < n_rph_sum_rep1; ++i){
    gamma_param[0][i] = std::pow(tmp1_sum[i], 2) * tph / tmp2_sum[i];
    gamma_param[1][i] = tmp2_sum[i] / tmp1_sum[i];
  }

  return(gamma_param);
}

// [[Rcpp::export]]
double gamma_pdf(double value, double shape, double scale){
  boost::math::gamma_distribution<double> d(shape, scale);
  return(pdf(d, value));
}

// [[Rcpp::export]]
double gamma_cdf(double threshold, double shape, double scale){
  boost::math::gamma_distribution<double> p(shape, scale);
  return(cdf(p, threshold));
}

// [[Rcpp::export]]
std::vector<double> calc_lpd_dedup1_all_mk(int n_dedup1_sum_rep1,
                                           int n_general_sum_rep1,
                                           int n_dropout_sum_rep1,
                                           int n_dropin_sum_rep1,
                                           std::vector<double> oph_rep1,
                                           std::vector<std::vector<double>> gamma_param,
                                           std::vector<double> at_dropout_rep1,
                                           std::vector<double> loglike_di){
  std::vector<double> lpd_dedup1_all_mk(n_dedup1_sum_rep1, 0);
  std::vector<double> gamma_shape =  gamma_param.at(0);
  std::vector<double> gamma_scale =  gamma_param.at(1);

  if(n_general_sum_rep1 != 0){
    std::vector<double> oph_general = extract_subvec_double(oph_rep1, 0, n_general_sum_rep1 - 1);
    std::vector<double> gamma_shape_general = extract_subvec_double(gamma_shape, 0, n_general_sum_rep1 - 1);
    std::vector<double> gamma_scale_general = extract_subvec_double(gamma_scale, 0, n_general_sum_rep1 - 1);
    for(int i = 0; i < n_general_sum_rep1; ++i){
      lpd_dedup1_all_mk[i] = log(gamma_pdf(oph_general[i], gamma_shape_general[i], gamma_scale_general[i]));
    }
  }

  if(n_dropout_sum_rep1 != 0){
    std::vector<double> gamma_shape_dropout = extract_subvec_double(gamma_shape, n_general_sum_rep1, n_general_sum_rep1 + n_dropout_sum_rep1 - 1);
    std::vector<double> gamma_scale_dropout = extract_subvec_double(gamma_scale, n_general_sum_rep1, n_general_sum_rep1 + n_dropout_sum_rep1 - 1);
    for(int i = 0; i < n_dropout_sum_rep1; ++i){
      lpd_dedup1_all_mk[n_general_sum_rep1 + i] = log(gamma_cdf(at_dropout_rep1[i], gamma_shape_dropout[i], gamma_scale_dropout[i]));
    }
  }

  if(n_dropin_sum_rep1 != 0){
    for(int i = 0; i < n_dropin_sum_rep1; ++i){
      lpd_dedup1_all_mk[n_general_sum_rep1 + n_dropout_sum_rep1 + i] = loglike_di[i];
    }
  }

  return(lpd_dedup1_all_mk);
}


std::vector<double> calc_lpd_gtc_cpp(int n_peak,
                                     int n_gtc,
                                     std::vector<double> lpd_dedup1_all_mk,
                                     std::vector<int> index_dedup2){
  std::vector<double> lpd_gtc(n_gtc, 0);
  for(int i = 0; i < n_gtc; ++i){
    for(int j = 0; j < n_peak; ++j){
      int index = index_dedup2[n_peak * i + j];
      lpd_gtc[i] += lpd_dedup1_all_mk[index];
    }
  }
  return(lpd_gtc);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_lpd_gtc_all_mk_cpp(std::vector<double> mp,
                                                         std::vector<double> deg,
                                                         int hnc,
                                                         int n_rph_sum_rep1,
                                                         std::vector<double> rph_base_use,
                                                         std::vector<std::vector<int>> index_dedup1_rep1_cpp,
                                                         std::vector<double> size_rep1,
                                                         double min_size,
                                                         std::vector<std::vector<double>> cv_square_rep1,
                                                         double tph,
                                                         int n_dedup1_sum_rep1,
                                                         int n_general_sum_rep1,
                                                         int n_dropout_sum_rep1,
                                                         int n_dropin_sum_rep1,
                                                         std::vector<double> oph_rep1,
                                                         std::vector<double> at_dropout_rep1,
                                                         std::vector<double> loglike_di,
                                                         int n_mk,
                                                         std::vector<int> index_dedup2_use_rep1_cpp,
                                                         std::vector<int> n_dedup2_use,
                                                         std::vector<int> n_peak_use,
                                                         std::vector<int> n_gtc_use){

  std::vector<std::vector<double>> rph_dedup1 = calc_rph_dedup1_cpp(mp, deg, hnc, n_rph_sum_rep1, rph_base_use, index_dedup1_rep1_cpp, size_rep1, min_size);

  std::vector<std::vector<double>> gamma_param = calc_gamma_param_cpp(rph_dedup1, cv_square_rep1, 2 * hnc, n_rph_sum_rep1, tph);

  std::vector<double> lpd_dedup1_all_mk = calc_lpd_dedup1_all_mk(n_dedup1_sum_rep1, n_general_sum_rep1, n_dropout_sum_rep1, n_dropin_sum_rep1,
                                                                 oph_rep1, gamma_param, at_dropout_rep1, loglike_di);

  std::vector<std::vector<double>> lpd_gtc_all_mk(n_mk);
  int index_start = 0;
  for(int i = 0; i < n_mk; ++i){
    int n_dedup2 = n_dedup2_use[i];
    int n_peak = n_peak_use[i];
    int n_gtc = n_gtc_use[i];
    std::vector<int> index_dedup2 = extract_subvec_int(index_dedup2_use_rep1_cpp, index_start, index_start + n_dedup2 - 1);
    lpd_gtc_all_mk.at(i) = calc_lpd_gtc_cpp(n_peak, n_gtc, lpd_dedup1_all_mk, index_dedup2);
    index_start += n_dedup2;
  }
  return(lpd_gtc_all_mk);
}
