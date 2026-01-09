#include "header.h"

// [[Rcpp::export]]
std::vector<double> calc_eph_dedup1_cpp(const std::vector<double>& temp_rep1,
                                        const std::vector<double>& deg,
                                        const int hnc,
                                        const std::vector<double>& rph_base_use_mk1,
                                        const std::vector<int>& index_dedup1,
                                        const double deg_size){
  std::vector<double> eph_dedup1(2 * hnc, 0);

  for(int i = 0; i < hnc; ++i){
    double temp_c1 = temp_rep1[i];
    double deg_c1 = deg[i];

    int index_al1 = index_dedup1[2 * i];
    int index_al2 = index_dedup1[2 * i + 1];

    double eph_base = temp_c1 * std::pow(deg_c1, deg_size);

    eph_dedup1[2 * i] = rph_base_use_mk1[index_al1] * eph_base;
    eph_dedup1[2 * i + 1] = rph_base_use_mk1[index_al2] * eph_base;
  }

  return(eph_dedup1);
}

// [[Rcpp::export]]
double gamma_pdf(const double value, const double shape, const double scale){

  // Check input values
  if(!std::isfinite(value) || value < 0.0 ||
     !std::isfinite(shape) || shape <= 0.0 ||
     !std::isfinite(scale) || scale <= 0.0){
    return 0.0;
  }

  try{
    boost::math::gamma_distribution<double> d(shape, scale);
    double result = pdf(d, value);
    if(!std::isfinite(result)){
      return 0.0;
    }
    return result;
  }catch(...){
    return 0.0;
  }
}

// [[Rcpp::export]]
double gamma_cdf(const double threshold, const double shape, const double scale){

  // Check input values
  if(!std::isfinite(threshold) || threshold < 0.0 ||
     !std::isfinite(shape) || shape <= 0.0 ||
     !std::isfinite(scale) || scale <= 0.0){
     return 0.0;
  }

  try{
    boost::math::gamma_distribution<double> p(shape, scale);
    double result = cdf(p, threshold);
    if(!std::isfinite(result)){
      return 0.0;
    }
    return result;
  }catch(...){
    return 0.0;
  }
}

// [[Rcpp::export]]
std::vector<double> calc_gamma_param_cpp(const std::vector<double>& eph_dedup1,
                                         const std::vector<double>& scale_base_mk1,
                                         const std::vector<int>& index_dedup1){

  double tmp1_sum(0);
  double tmp2_sum(0);
  int n_al = eph_dedup1.size();

  for(int i = 0; i < n_al; ++i){
    tmp1_sum += eph_dedup1[i];
    tmp2_sum += eph_dedup1[i] * scale_base_mk1[index_dedup1[i]];
  }

  std::vector<double> gamma_param(2);
  gamma_param[0] = std::pow(tmp1_sum, 2) / tmp2_sum;
  gamma_param[1] = tmp2_sum / tmp1_sum;

  return(gamma_param);
}

// [[Rcpp::export]]
double calc_gamma_cdf_dropout(const std::vector<double>& eph_dedup1,
                              const std::vector<double>& scale_base_mk1,
                              const std::vector<int>& index_dedup1,
                              const int index_ignore,
                              const double at){

  int n_al = eph_dedup1.size();
  double log_cd_sum = 0;

  for(int i = 0; i < n_al; ++i){
    int index_target = index_dedup1[i];

    if(index_target != index_ignore){
      double gamma_shape = eph_dedup1[i] / scale_base_mk1[index_target];
      double gamma_scale = scale_base_mk1[index_target];
      double cd = gamma_cdf(at, gamma_shape, gamma_scale);

      double log_cd;
      if(cd > 0.0 && std::isfinite(cd)){
        log_cd = std::log(cd);
      }else{
        log_cd = -std::numeric_limits<double>::infinity();
      }
      log_cd_sum += log_cd;
    }
  }

  return log_cd_sum;
}

// [[Rcpp::export]]
std::vector<double> calc_lpd_dedup1_mk(const std::vector<double>& temp_rep1,
                                       const std::vector<double>& deg,
                                       const int hnc,
                                       const double at,
                                       const int n_general,
                                       const int n_dropout,
                                       const int n_dropin,
                                       const int n_ignore,
                                       const std::vector<double>& rph_base_use_mk1,
                                       const std::vector<double>& scale_base_mk1,
                                       const std::vector<double>& deg_size_base_mk1,
                                       const std::vector<double>& oph_base_rep1_mk1,
                                       const std::vector<double>& loglike_di_rep1_mk1,
                                       const std::vector<std::vector<int>>& index_dedup1_rep1_mk1_cpp){

  std::vector<double> lpd_dedup1_mk(n_general + n_dropout + n_dropin + n_ignore, 0);

  if(n_general != 0){
    #pragma omp parallel for
    for(int i = 0; i < n_general; ++i){
      std::vector<int> index_dedup1 = index_dedup1_rep1_mk1_cpp[i];
      int min_index = *std::min_element(index_dedup1.begin(), index_dedup1.end());

      std::vector<double> eph_dedup1_general = calc_eph_dedup1_cpp(temp_rep1, deg, hnc, rph_base_use_mk1, index_dedup1, deg_size_base_mk1[min_index]);
      std::vector<double> gamma_param = calc_gamma_param_cpp(eph_dedup1_general, scale_base_mk1, index_dedup1);

      double pd = gamma_pdf(oph_base_rep1_mk1[min_index], gamma_param[0], gamma_param[1]);
      double log_pd;
      if(pd > 0.0 && std::isfinite(pd)){
        log_pd = std::log(pd);
      }else{
        log_pd = -std::numeric_limits<double>::infinity();
      }
      lpd_dedup1_mk[i] = log_pd;
    }
  }

  if(n_dropout != 0){
    int index_ignore = rph_base_use_mk1.size() - 2;

    #pragma omp parallel for
    for(int i = 0; i < n_dropout; ++i){
      std::vector<int> index_dedup1 = index_dedup1_rep1_mk1_cpp[n_general + i];
      int min_index = *std::min_element(index_dedup1.begin(), index_dedup1.end());

      std::vector<double> eph_dedup1_dropout = calc_eph_dedup1_cpp(temp_rep1, deg, hnc, rph_base_use_mk1, index_dedup1, deg_size_base_mk1[min_index]);

      lpd_dedup1_mk[n_general + i] = calc_gamma_cdf_dropout(eph_dedup1_dropout, scale_base_mk1, index_dedup1, index_ignore, at);
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
std::vector<std::vector<double>> calc_lpd_gtc_rep1(const std::vector<double>& temp_rep1,
                                                   const std::vector<double>& deg,
                                                   const int hnc,
                                                   const std::vector<double>& ats_use,
                                                   const std::vector<int>& n_general_rep1,
                                                   const std::vector<int>& n_dropout_rep1,
                                                   const std::vector<int>& n_dropin_rep1,
                                                   const std::vector<int>& n_ignore_rep1,
                                                   const std::vector<std::vector<double>>& rph_base_use,
                                                   const std::vector<std::vector<double>>& scale_base_use,
                                                   const std::vector<std::vector<double>>& deg_size_base_use,
                                                   const std::vector<std::vector<double>>& oph_base_use_rep1,
                                                   const std::vector<std::vector<double>>& loglike_di_rep1,
                                                   const std::vector<std::vector<std::vector<int>>>& index_dedup1_rep1_cpp,
                                                   const std::vector<std::vector<int>>& index_dedup2_rep1_cpp,
                                                   const std::vector<int>& n_gtc_use,
                                                   const std::vector<int>& n_peak_use){

  int n_mk = n_general_rep1.size();
  std::vector<std::vector<double>> lpd_gtc_rep1(n_mk, std::vector<double>());

  for(int i = 0; i < n_mk; ++i){
    std::vector<double> lpd_dedup1_mk = calc_lpd_dedup1_mk(temp_rep1, deg, hnc, ats_use[i], n_general_rep1[i], n_dropout_rep1[i], n_dropin_rep1[i], n_ignore_rep1[i],
                                                           rph_base_use[i], scale_base_use[i], deg_size_base_use[i], oph_base_use_rep1[i], loglike_di_rep1[i],
                                                           index_dedup1_rep1_cpp[i]);

    lpd_gtc_rep1[i] = calc_lpd_gtc_mk(lpd_dedup1_mk, index_dedup2_rep1_cpp[i], n_gtc_use[i], n_peak_use[i]);
  }

  return(lpd_gtc_rep1);
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> calc_lpd_gtc_all_rep(const std::vector<double>& temp,
                                                                   const std::vector<double>& deg,
                                                                   const std::vector<double>& rep_ae,
                                                                   const std::vector<double>& ats_use,
                                                                   const std::vector<std::vector<int>>& n_general_list,
                                                                   const std::vector<std::vector<int>>& n_dropout_list,
                                                                   const std::vector<std::vector<int>>& n_dropin_list,
                                                                   const std::vector<std::vector<int>>& n_ignore_list,
                                                                   const std::vector<std::vector<double>>& rph_base_use,
                                                                   const std::vector<std::vector<double>>& scale_base_use,
                                                                   const std::vector<std::vector<double>>& deg_size_base_use,
                                                                   const std::vector<std::vector<std::vector<double>>>& oph_base_use,
                                                                   const std::vector<std::vector<std::vector<double>>>& loglike_di_list,
                                                                   const std::vector<std::vector<std::vector<std::vector<int>>>>& index_dedup1_list_cpp,
                                                                   const std::vector<std::vector<std::vector<int>>>& index_dedup2_list_cpp,
                                                                   const std::vector<int>& n_gtc_use,
                                                                   const std::vector<int>& n_peak_use){

  int n_rep = rep_ae.size();
  int hnc = temp.size();
  int n_mk = n_gtc_use.size();
  std::vector<std::vector<std::vector<double>>> lpd_gtc_all_rep(n_rep, std::vector<std::vector<double>>(n_mk, std::vector<double>()));

  for(int i = 0; i < n_rep; ++i){

    std::vector<double> temp_rep1 = temp;
    for(int j = 0; j < hnc; ++j){
      temp_rep1[j] *= rep_ae[i];
    }

    lpd_gtc_all_rep[i] = calc_lpd_gtc_rep1(temp_rep1, deg, hnc, ats_use, n_general_list[i], n_dropout_list[i], n_dropin_list[i], n_ignore_list[i],
                                           rph_base_use, scale_base_use, deg_size_base_use, oph_base_use[i], loglike_di_list[i],
                                           index_dedup1_list_cpp[i], index_dedup2_list_cpp[i], n_gtc_use, n_peak_use);
  }

  return(lpd_gtc_all_rep);
}

// [[Rcpp::export]]
Rcpp::List calc_loglike_cpp(const std::vector<std::vector<std::vector<double>>>& lpd_gtc_all_rep,
                              const std::vector<std::vector<double>>& gtc_logprob_use){

  int n_mk = gtc_logprob_use.size();
  int n_rep = lpd_gtc_all_rep.size();
  std::vector<double> loglike_all_mk(n_mk, 0);
  std::vector<std::vector<double>> loglike_gtc_all_mk(n_mk, std::vector<double>());

  for(int i = 0; i < n_mk; ++i){
    std::vector<double> loglike_gtc_mk = gtc_logprob_use[i];
    int n_gtc = loglike_gtc_mk.size();

    for(int j = 0; j < n_rep; ++j){
      std::vector<double> lpd_gtc_mk = lpd_gtc_all_rep[j][i];

      for(int k = 0; k < n_gtc; ++k){
        loglike_gtc_mk[k] += lpd_gtc_mk[k];
      }
    }

    loglike_all_mk[i] = log_sum_exp(loglike_gtc_mk);
    loglike_gtc_all_mk[i] = loglike_gtc_mk;
  }

  List result = List::create(loglike_all_mk, loglike_gtc_all_mk);
  return(result);
}
