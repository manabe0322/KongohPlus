#include "header.h"

// [[Rcpp::export]]
double gamma_pdf(double value, double alpha, double beta){
  boost::math::gamma_distribution<double> d(alpha, 1 / beta);
  return(pdf(d, value));
}

// [[Rcpp::export]]
double gamma_cdf(double threshold, double alpha, double beta){
  boost::math::gamma_distribution<double> p(alpha, 1 / beta);
  return(cdf(p, threshold));
}

// [[Rcpp::export]]
double calc_lpd(double tph, double oph, double eph, double at, double beta, double di_alpha, double di_beta, double logprob_di, double di_min){
  double lpd = 0;
  if(eph > 0){
    if(oph >= at){
      lpd = log(gamma_pdf(oph, std::pow(beta, -2) * eph, 1 / (tph * std::pow(beta, 2))));
    }else{
      lpd = log(gamma_cdf(at, std::pow(beta, -2) * eph, 1 / (tph * std::pow(beta, 2))));
    }
  //drop-in
  }else if(eph < 0){
    lpd = log(gamma_pdf(oph - di_min, di_alpha, di_beta)) + logprob_di;
  }
  return(lpd);
}

// [[Rcpp::export]]
std::vector<double> calc_lpd_mk(double tph,
                                std::vector<double> oph_mk,
                                std::vector<double> eph_mk,
                                double at,
                                double beta,
                                double di_alpha,
                                double di_beta,
                                double logprob_di,
                                double di_min){
  int n_eph = eph_mk.size();
  std::vector<double> lpd_mk(n_eph);
  for(int i = 0; i < n_eph; ++i){
    lpd_mk[i] = calc_lpd(tph, oph_mk[i], eph_mk[i], at, beta, di_alpha, di_beta, logprob_di, di_min);
  }
  return(lpd_mk);
}

// [[Rcpp::export]]
std::vector<double> calc_lpd_gtc(double tph,
                                 std::vector<double> oph_mk,
                                 std::vector<double> eph_mk,
                                 double at,
                                 double beta,
                                 double di_alpha,
                                 double di_beta,
                                 double logprob_di,
                                 double di_min,
                                 int n_peak,
                                 int n_gtc,
                                 std::vector<int> index_dedup2_mk){
  std::vector<double> lpd_mk = calc_lpd_mk(tph, oph_mk, eph_mk, at, beta, di_alpha, di_beta, logprob_di, di_min);
  std::vector<double> lpd_gtc(n_gtc, 0);

  int index_start = 0;
  for(int i = 0; i < n_gtc; ++i){
    std::vector<int> index_dedup2_gtc = extract_subvec_int(index_dedup2_mk, index_start, index_start + n_peak - 1);
    for(int j = 0; j < n_peak; ++j){
      lpd_gtc[i] += lpd_mk[index_dedup2_gtc[j]];
    }
    index_start += n_peak;
  }
  return(lpd_gtc);
}

// [[Rcpp::export]]
std::vector<double> calc_eph_mk(std::vector<double> mp,
                                std::vector<double> eph_base_mk,
                                std::vector<std::vector<int>> index_dedup1_mk){
  int hnc = mp.size();
  int n_dedup1 = index_dedup1_mk.at(0).size();
  std::vector<double> eph_mk(n_dedup1, 0);
  for(int i = 0; i < hnc; ++i){
    double mp_c1 = mp[i];
    std::vector<int> ind_dedup1_al1 = index_dedup1_mk.at(2 * i);
    std::vector<int> ind_dedup1_al2 = index_dedup1_mk.at(2 * i + 1);
    for(int j = 0; j < n_dedup1; ++j){
      eph_mk[j] += eph_base_mk[ind_dedup1_al1[j]] * mp_c1;
      eph_mk[j] += eph_base_mk[ind_dedup1_al2[j]] * mp_c1;
    }
  }
  return(eph_mk);
}

std::vector<double> calc_loglike_gtc(int n_gtc,
                                     std::vector<double> lpd_gtc,
                                     std::vector<double> gtc_logprob_mk){
  std::vector<double> loglike_gtc = gtc_logprob_mk;
  for(int i = 0; i < n_gtc; ++i){
    loglike_gtc[i] += lpd_gtc[i];
  }
  return(loglike_gtc);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_loglike_gtc_all_mp(double oph_sum,
                                                         std::vector<std::vector<double>> mp_cand,
                                                         std::vector<double> eph_base_mk,
                                                         std::vector<std::vector<int>> index_dedup1_mk,
                                                         std::vector<double> oph_mk,
                                                         double at,
                                                         double beta,
                                                         double di_alpha,
                                                         double di_beta,
                                                         double logprob_di,
                                                         double di_min,
                                                         int n_peak,
                                                         int n_gtc,
                                                         std::vector<int> index_dedup2_mk,
                                                         std::vector<double> gtc_logprob_mk){
  int n_mp_cand = mp_cand.size();
  std::vector<std::vector<double>> loglike_gtc_all_mp(n_mp_cand, std::vector<double>(n_gtc));
  for(int i = 0; i < n_mp_cand; ++i){
    std::vector<double> mp = mp_cand.at(i);
    std::vector<double> eph_mk = calc_eph_mk(mp, eph_base_mk, index_dedup1_mk);
    std::vector<double> lpd_gtc = calc_lpd_gtc(oph_sum, oph_mk, eph_mk, at, beta, di_alpha, di_beta, logprob_di, di_min, n_peak, n_gtc, index_dedup2_mk);
    std::vector<double> loglike_gtc = calc_loglike_gtc(n_gtc, lpd_gtc, gtc_logprob_mk);
    loglike_gtc_all_mp.at(i) = loglike_gtc;
  }
  return(loglike_gtc_all_mp);
}

// [[Rcpp::export]]
std::vector<double> calc_eph_dedup1_cpp(std::vector<double> mp,
                                        std::vector<double> deg,
                                        int hnc,
                                        int n_dedup1,
                                        std::vector<double> eph_base,
                                        std::vector<std::vector<int>> index_dedup1,
                                        std::vector<double> size_dedup1){
  std::vector<double> eph_dedup1(n_dedup1, 0);
  for(int i = 0; i < hnc; ++i){
    double mp_c1 = mp[i];
    double deg_c1 = deg[i];
    std::vector<int> ind_dedup1_al1 = index_dedup1.at(2 * i);
    std::vector<int> ind_dedup1_al2 = index_dedup1.at(2 * i + 1);
    for(int j = 0; j < n_dedup1; ++j){
      eph_dedup1[j] += eph_base[ind_dedup1_al1[j]] * mp_c1 * exp(deg_c1 * size_dedup1[j]);
      eph_dedup1[j] += eph_base[ind_dedup1_al2[j]] * mp_c1 * exp(deg_c1 * size_dedup1[j]);
    }
  }
  return(eph_dedup1);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_lpd_gtc_adopt(int n_mk,
                                                    double tph,
                                                    std::vector<int> n_dedup1_adopt,
                                                    std::vector<double> oph_dedup1,
                                                    std::vector<double> eph_dedup1,
                                                    std::vector<double> at_mk,
                                                    double beta,
                                                    double di_alpha,
                                                    double di_beta,
                                                    double logprob_di,
                                                    double di_min,
                                                    std::vector<int> n_peak_adopt,
                                                    std::vector<int> n_gtc_adopt,
                                                    std::vector<std::vector<int>> index_dedup2_adopt_cpp){
  std::vector<std::vector<double>> lpd_gtc_adopt(n_mk, std::vector<double>());
  int eph_start = 0;
  for(int i = 0; i < n_mk; ++i){
    int n_dedup_eph = n_dedup1_adopt[i];
    std::vector<double> oph_mk = extract_subvec_double(oph_dedup1, eph_start, eph_start + n_dedup_eph - 1);
    std::vector<double> eph_mk = extract_subvec_double(eph_dedup1, eph_start, eph_start + n_dedup_eph - 1);
    double at = at_mk[i];
    int n_peak = n_peak_adopt[i];
    int n_gtc = n_gtc_adopt[i];

    std::vector<int> index_dedup2_mk = index_dedup2_adopt_cpp.at(i);
    lpd_gtc_adopt.at(i) = calc_lpd_gtc(tph, oph_mk, eph_mk, at, beta, di_alpha, di_beta, logprob_di, di_min, n_peak, n_gtc, index_dedup2_mk);

    eph_start += n_dedup_eph;
  }
  return(lpd_gtc_adopt);
}
