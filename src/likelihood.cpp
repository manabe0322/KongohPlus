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

std::vector<double> calc_lpd_general_cpp(int n_general_sum,
                                         std::vector<double> oph_general,
                                         std::vector<double> rph_general,
                                         double cv_m2,
                                         double cv_p2,
                                         double tph){
  std::vector<double> lpd_general(n_general_sum);
  double gamma_scale = tph * cv_p2;
  for(int i = 0; i < n_general_sum, ++i){
    lpd_general[i] = log(gamma_pdf(oph_general[i], cv_m2 * rph_general[i], gamma_scale));
  }
  return(lpd_general);
}

std::vector<double> calc_lpd_dropout_cpp(int n_dropout_sum,
                                         std::vector<double> at_dropout,
                                         std::vector<double> rph_dropout,
                                         double cv_m2,
                                         double cv_p2,
                                         double tph){
  std::vector<double> lpd_dropout(n_dropout_sum);
  double gamma_scale = tph * cv_p2;
  for(int i = 0; i < n_dropout_sum, ++i){
    lpd_dropout[i] = log(gamma_cdf(at_dropout[i], cv_m2 * rph_dropout[i], gamma_scale));
  }
  return(lpd_dropout);
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
    for(int i = 0; i < n_general_sum, ++i){
      lpd_dedup1_all_mk[i] = log(gamma_pdf(oph_general[i], cv_m2 * rph_general[i], gamma_scale));
    }
  }

  if(n_dropout_sum != 0){
    std::vector<double> rph_dropout = extract_subvec_double(rph_dedup1_all_mk, n_general_sum, n_general_sum + n_dropout_sum - 1);
    for(int i = 0; i < n_dropout_sum, ++i){
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
                                     std::vector<double> lpd_dedup2,
                                     std::vector<double> gtc_logprob_mk){
  std::vector<double> lpd_gtc(n_gtc, 0);
  for(int i = 0; i < n_gtc; ++i){
    std::vector<double> lpd_dedup2_gtc1 = extract_subvec_double(lpd_dedup2, n_peak * i, n_peak * (i + 1) - 1);
    lpd_gtc[i] += accumulate(lpd_dedup2_gtc1.begin(), lpd_dedup2_gtc1.end(), 0.0);
  }
  return(lpd_gtc);
}

std::vector<double> calc_loglike_gtc_cpp(int n_gtc,
                                         std::vector<double> lpd_gtc,
                                         std::vector<double> gtc_logprob_mk){
  std::vector<double> loglike_gtc = gtc_logprob_mk;
  for(int i = 0; i < n_gtc; ++i){
    loglike_gtc[i] += lpd_gtc[i];
  }
  return(loglike_gtc);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_loglike_cpp(int n_mk,
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

  int index_start = 0;
  for(int i = 0; i < n_mk; ++i){
    int n_dedup2 = n_dedup2_all_mk[i];
    int n_peak = n_peak_all_mk[i];
    int n_gtc = n_gtc_all_mk[i];
    std::vector<int> index_dedup2 = extract_subvec_int(index_dedup2_all_mk_cpp, index_start, index_start + n_dedup2 - 1);

    index_start += n_dedup2;
  }
}
