#include "header.h"

// [[Rcpp::export]]
std::vector<double> calc_rph_dedup1_cpp(std::vector<double> mp,
                                        std::vector<double> deg,
                                        int hnc,
                                        int n_dedup1_sum,
                                        std::vector<double> rph_base_all_mk,
                                        std::vector<std::vector<int>> index_dedup1_all_mk_cpp,
                                        std::vector<double> size_dedup1_all_mk,
                                        std::vector<double> ae_dedup1_all_mk){
  std::vector<double> rph_dedup1_all_mk(n_dedup1_sum, 0);
  for(int i = 0; i < hnc; ++i){
    double mp_c1 = mp[i];
    double deg_c1 = deg[i];
    std::vector<int> ind_dedup1_al1 = index_dedup1_all_mk_cpp.at(2 * i);
    std::vector<int> ind_dedup1_al2 = index_dedup1_all_mk_cpp.at(2 * i + 1);
    for(int j = 0; j < n_dedup1_sum; ++j){
      rph_dedup1_all_mk[j] += rph_base_all_mk[ind_dedup1_al1[j]] * mp_c1 * exp(deg_c1 * size_dedup1_all_mk[j]) * ae_dedup1_all_mk[j];
      rph_dedup1_all_mk[j] += rph_base_all_mk[ind_dedup1_al2[j]] * mp_c1 * exp(deg_c1 * size_dedup1_all_mk[j]) * ae_dedup1_all_mk[j];
    }
  }
  return(rph_dedup1_all_mk);
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

std::vector<double> calc_lpd_dedup1_all_mk(int n_dedup1_sum,
                                           int n_general_sum,
                                           int n_dropout_sum,
                                           int n_dropin_sum,
                                           std::vector<double> oph_dedup1_all_mk,
                                           std::vector<double> rph_dedup1_all_mk,
                                           std::vector<double> at_dropout,
                                           std::vector<double> loglike_di,
                                           double cv,
                                           double tph){
  std::vector<double> lpd_dedup1_all_mk(n_dedup1_sum, 0);
  double cv_m2 = std::pow(cv, -2);
  double gamma_scale = tph * std::pow(cv, 2);

  if(n_general_sum != 0){
    std::vector<double> oph_general = extract_subvec_double(oph_dedup1_all_mk, 0, n_general_sum - 1);
    std::vector<double> rph_general = extract_subvec_double(rph_dedup1_all_mk, 0, n_general_sum - 1);
    for(int i = 0; i < n_general_sum; ++i){
      lpd_dedup1_all_mk[i] = log(gamma_pdf(oph_general[i], cv_m2 * rph_general[i], gamma_scale));
    }
  }

  if(n_dropout_sum != 0){
    std::vector<double> rph_dropout = extract_subvec_double(rph_dedup1_all_mk, n_general_sum, n_general_sum + n_dropout_sum - 1);
    for(int i = 0; i < n_dropout_sum; ++i){
      lpd_dedup1_all_mk[n_general_sum + i] = log(gamma_cdf(at_dropout[i], cv_m2 * rph_dropout[i], gamma_scale));
    }
  }

  if(n_dropin_sum != 0){
    for(int i = 0; i < n_dropin_sum; ++i){
      lpd_dedup1_all_mk[n_general_sum + n_dropout_sum + i] = loglike_di[i];
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
std::vector<std::vector<double>> calc_lpd_gtc_all_mk_cpp(int n_mk,
                                                         std::vector<double> mp,
                                                         std::vector<double> deg,
                                                         int hnc,
                                                         int n_dedup1_sum,
                                                         int n_general_sum,
                                                         int n_dropout_sum,
                                                         int n_dropin_sum,
                                                         std::vector<double> rph_base_all_mk,
                                                         std::vector<std::vector<int>> index_dedup1_all_mk_cpp,
                                                         std::vector<double> size_dedup1_all_mk,
                                                         std::vector<double> ae_dedup1_all_mk,
                                                         std::vector<double> oph_dedup1_all_mk,
                                                         std::vector<double> at_dropout,
                                                         std::vector<double> loglike_di,
                                                         double cv,
                                                         double tph,
                                                         std::vector<int> index_dedup2_all_mk_cpp,
                                                         std::vector<int> n_dedup2_all_mk,
                                                         std::vector<int> n_peak_all_mk,
                                                         std::vector<int> n_gtc_all_mk){

  std::vector<double> rph_dedup1_all_mk = calc_rph_dedup1_cpp(mp, deg, hnc, n_dedup1_sum, rph_base_all_mk, index_dedup1_all_mk_cpp, size_dedup1_all_mk, ae_dedup1_all_mk);

  std::vector<double> lpd_dedup1_all_mk = calc_lpd_dedup1_all_mk(n_dedup1_sum, n_general_sum, n_dropout_sum, n_dropin_sum,
                                                                 oph_dedup1_all_mk, rph_dedup1_all_mk, at_dropout, loglike_di, cv, tph);

  std::vector<std::vector<double>> lpd_gtc_all_mk(n_mk);
  int index_start = 0;
  for(int i = 0; i < n_mk; ++i){
    int n_dedup2 = n_dedup2_all_mk[i];
    int n_peak = n_peak_all_mk[i];
    int n_gtc = n_gtc_all_mk[i];
    std::vector<int> index_dedup2 = extract_subvec_int(index_dedup2_all_mk_cpp, index_start, index_start + n_dedup2 - 1);
    lpd_gtc_all_mk.at(i) = calc_lpd_gtc_cpp(n_peak, n_gtc, lpd_dedup1_all_mk, index_dedup2);
    index_start += n_dedup2;
  }
  return(lpd_gtc_all_mk);
}
