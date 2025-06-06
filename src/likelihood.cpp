#include "header.h"

// [[Rcpp::export]]
std::vector<double> calc_rph_dedup1_cpp(std::vector<double> mp,
                                        std::vector<double> deg,
                                        int hnc,
                                        std::vector<double> rph_base_use,
                                        std::vector<int> index_dedup1,
                                        double deg_size){
  std::vector<double> rph_dedup1(2 * hnc, 0);

  for(int i = 0; i < hnc; ++i){
    double mp_c1 = mp[i];
    double deg_c1 = deg[i];

    int index_al1 = index_dedup1[2 * i];
    int index_al2 = index_dedup1[2 * i + 1];

    double rph_mpdeg = mp_c1 * std::pow(deg_c1, deg_size);

    rph_dedup1[2 * i] = rph_base_use[index_al1] * rph_mpdeg;
    rph_dedup1[2 * i + 1] = rph_base_use[index_al2] * rph_mpdeg;
  }

  return(rph_dedup1);
}

// [[Rcpp::export]]
std::vector<double> calc_gamma_param_cpp(std::vector<double> rph_dedup1,
                                         std::vector<double> cv_squared_dedup1,
                                         double tph){

  double tmp1_sum(0);
  double tmp2_sum(0);
  int n_al = rph_dedup1.size();

  for(int i = 0; i < n_al; ++i){
    tmp1_sum += rph_dedup1[i];
    tmp2_sum += rph_dedup1[i] * cv_squared_dedup1[i];
  }

  std::vector<double> gamma_param(2);
  gamma_param[0] = std::pow(tmp1_sum, 2) * tph / tmp2_sum;
  gamma_param[1] = tmp2_sum / tmp1_sum;

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
std::vector<double> calc_lpd_dedup1_mk(double tph,
                                       std::vector<double> mp,
                                       std::vector<double> deg,
                                       int hnc,
                                       double at,
                                       int n_general,
                                       int n_dropout,
                                       int n_dropin,
                                       int n_ignore,
                                       std::vector<double> rph_base_use_mk1,
                                       std::vector<std::vector<int>> index_dedup1_rep1_mk1_cpp,
                                       std::vector<std::vector<double>> cv_squared_dedup1_rep1_mk1,
                                       std::vector<double> deg_size_dedup1_rep1_mk1,
                                       std::vector<double> oph_dedup1_rep1_mk1,
                                       std::vector<double> loglike_di_rep1_mk1){

  std::vector<double> lpd_dedup1_mk(n_general + n_dropout + n_dropin + n_ignore, 0);

  if(n_general != 0){
    #pragma omp parallel for
    for(int i = 0; i < n_general; ++i){
      std::vector<double> rph_dedup1_general = calc_rph_dedup1_cpp(mp, deg, hnc, rph_base_use_mk1, index_dedup1_rep1_mk1_cpp[i], deg_size_dedup1_rep1_mk1[i]);
      std::vector<double> gamma_param = calc_gamma_param_cpp(rph_dedup1_general, cv_squared_dedup1_rep1_mk1[i], tph);
      lpd_dedup1_mk[i] = log(gamma_pdf(oph_dedup1_rep1_mk1[i], gamma_param[0], gamma_param[1]));
    }
  }

  if(n_dropout != 0){
    #pragma omp parallel for
    for(int i = 0; i < n_dropout; ++i){
      std::vector<double> rph_dedup1_dropout = calc_rph_dedup1_cpp(mp, deg, hnc, rph_base_use_mk1, index_dedup1_rep1_mk1_cpp[n_general + i], deg_size_dedup1_rep1_mk1[n_general + i]);
      std::vector<double> gamma_param = calc_gamma_param_cpp(rph_dedup1_dropout, cv_squared_dedup1_rep1_mk1[n_general + i], tph);
      lpd_dedup1_mk[n_general + i] = log(gamma_cdf(at, gamma_param[0], gamma_param[1]));
    }
  }

  if(n_dropin != 0){
    for(int i = 0; i < n_dropin; ++i){
      lpd_dedup1_mk[n_general + n_dropout + i] = loglike_di_rep1_mk1[i];
    }
  }

  return(lpd_dedup1_mk);
}

// [[Rcpp::export]]
std::vector<double> calc_lpd_gtc_mk(const std::vector<double>& lpd_dedup1_mk,
                                    const std::vector<int>& index_dedup2_rep1_mk1_cpp,
                                    const int n_gtc,
                                    const int n_peak){

  std::vector<double> lpd_gtc_mk(n_gtc, 0);

  #pragma omp parallel for
  for(int i = 0; i < n_gtc; ++i){
    for(int j = 0; j < n_peak; ++j){
      int index = index_dedup2_rep1_mk1_cpp[n_peak * i + j];
      lpd_gtc_mk[i] += lpd_dedup1_mk[index];
    }
  }

  return(lpd_gtc_mk);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_lpd_gtc_rep1(double tph,
                                                   std::vector<double> mp,
                                                   std::vector<double> deg,
                                                   int hnc,
                                                   std::vector<double> ats_use,
                                                   std::vector<int> n_general_rep1,
                                                   std::vector<int> n_dropout_rep1,
                                                   std::vector<int> n_dropin_rep1,
                                                   std::vector<int> n_ignore_rep1,
                                                   std::vector<std::vector<double>> rph_base_use,
                                                   std::vector<std::vector<std::vector<int>>> index_dedup1_rep1_cpp,
                                                   std::vector<std::vector<std::vector<double>>> cv_squared_dedup1_rep1,
                                                   std::vector<std::vector<double>> deg_size_dedup1_rep1,
                                                   std::vector<std::vector<double>> oph_dedup1_rep1,
                                                   std::vector<std::vector<double>> loglike_di_rep1,
                                                   std::vector<std::vector<int>> index_dedup2_rep1_cpp,
                                                   std::vector<int> n_gtc_use,
                                                   std::vector<int> n_peak_use){

  int n_mk = n_general_rep1.size();
  std::vector<std::vector<double>> lpd_gtc_rep1(n_mk, std::vector<double>());

  for(int i = 0; i < n_mk; ++i){
    std::vector<double> lpd_dedup1_mk = calc_lpd_dedup1_mk(tph, mp, deg, hnc, ats_use[i], n_general_rep1[i], n_dropout_rep1[i], n_dropin_rep1[i], n_ignore_rep1[i], rph_base_use[i],
                                                           index_dedup1_rep1_cpp[i], cv_squared_dedup1_rep1[i], deg_size_dedup1_rep1[i], oph_dedup1_rep1[i], loglike_di_rep1[i]);

    lpd_gtc_rep1[i] = calc_lpd_gtc_mk(lpd_dedup1_mk, index_dedup2_rep1_cpp[i], n_gtc_use[i], n_peak_use[i]);
  }

  return(lpd_gtc_rep1);
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> calc_lpd_gtc_all_rep(double tph,
                                                                   std::vector<double> mp,
                                                                   std::vector<double> deg,
                                                                   int hnc,
                                                                   std::vector<double> ats_use,
                                                                   std::vector<std::vector<int>> n_general_list,
                                                                   std::vector<std::vector<int>> n_dropout_list,
                                                                   std::vector<std::vector<int>> n_dropin_list,
                                                                   std::vector<std::vector<int>> n_ignore_list,
                                                                   std::vector<std::vector<double>> rph_base_use,
                                                                   std::vector<std::vector<std::vector<std::vector<int>>>> index_dedup1_list_cpp,
                                                                   std::vector<std::vector<std::vector<std::vector<double>>>> cv_squared_dedup1_list,
                                                                   std::vector<std::vector<std::vector<double>>> deg_size_dedup1_list,
                                                                   std::vector<std::vector<std::vector<double>>> oph_dedup1_list,
                                                                   std::vector<std::vector<std::vector<double>>> loglike_di_list,
                                                                   std::vector<std::vector<std::vector<int>>> index_dedup2_list_cpp,
                                                                   std::vector<int> n_gtc_use,
                                                                   std::vector<int> n_peak_use){

  int n_rep = n_general_list.size();
  int n_mk = n_gtc_use.size();
  std::vector<std::vector<std::vector<double>>> lpd_gtc_all_rep(n_rep, std::vector<std::vector<double>>(n_mk, std::vector<double>()));

  for(int i = 0; i < n_rep; ++i){
    lpd_gtc_all_rep[i] = calc_lpd_gtc_rep1(tph, mp, deg, hnc, ats_use, n_general_list[i], n_dropout_list[i], n_dropin_list[i], n_ignore_list[i], rph_base_use,
                                           index_dedup1_list_cpp[i], cv_squared_dedup1_list[i], deg_size_dedup1_list[i], oph_dedup1_list[i], loglike_di_list[i],
                                           index_dedup2_list_cpp[i], n_gtc_use, n_peak_use);
  }

  return(lpd_gtc_all_rep);
}

// [[Rcpp::export]]
std::vector<double> calc_loglike_cpp(const std::vector<std::vector<std::vector<double>>>& lpd_gtc_all_rep,
                                     const std::vector<std::vector<double>>& gtc_logprob_use){

  int n_mk = gtc_logprob_use.size();
  int n_rep = lpd_gtc_all_rep.size();
  std::vector<double> loglike_cpp(n_mk, 0);

  for(int i = 0; i < n_mk; ++i){
    std::vector<double> loglike_cpp_mk = gtc_logprob_use[i];
    int n_gtc = loglike_cpp_mk.size();

    for(int j = 0; j < n_rep; ++j){
      std::vector<double> lpd_gtc_mk = lpd_gtc_all_rep[j][i];

      for(int k = 0; k < n_gtc; ++k){
        loglike_cpp_mk[k] += lpd_gtc_mk[k];
      }
    }

    loglike_cpp[i] = log_sum_exp(loglike_cpp_mk);
  }

  return(loglike_cpp);
}
