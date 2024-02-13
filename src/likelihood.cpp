#include "header.h"

// [[Rcpp::export]]
std::vector<double> calc_eph_dedup_1_cpp(double tph,
                                         std::vector<double> mp,
                                         std::vector<double> deg,
                                         int hnc,
                                         int n_dedup_1,
                                         std::vector<double> eph_base,
                                         std::vector<std::vector<int>> index_dedup_1,
                                         std::vector<double> size_dedup_1){
  std::vector<double> eph_dedup_1(n_dedup_1, 0);
  for(int i = 0; i < hnc; ++i){
    double tph_c1 = tph * mp[i];
    double deg_c1 = deg[i];
    std::vector<int> ind_dedup1_al1 = index_dedup_1.at(2 * i);
    std::vector<int> ind_dedup1_al2 = index_dedup_1.at(2 * i + 1);
    for(int j = 0; j < n_dedup_1; ++j){
      eph_dedup_1[j] += eph_base[ind_dedup1_al1[j]] * tph_c1 * exp(deg_c1 * size_dedup_1[j]);
      eph_dedup_1[j] += eph_base[ind_dedup1_al2[j]] * tph_c1 * exp(deg_c1 * size_dedup_1[j]);
    }
  }
  return(eph_dedup_1);
}

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
double calc_lpd(double oph, double eph, int at, double beta){
  double lpd = 0;
  double alpha = eph * beta;
  if(eph != 0){
    if(oph >= at){
      lpd = log(gamma_pdf(oph, alpha, beta));
    }else{
      lpd = log(gamma_cdf(at, alpha, beta));
    }
  }
  return(lpd);
}

// [[Rcpp::export]]
std::vector<double> calc_lpd_mk(std::vector<double> oph_mk,
                                std::vector<double> eph_mk,
                                int at,
                                double beta){
  int n_eph = eph_mk.size();
  std::vector<double> lpd_mk(n_eph);
  for(int i = 0; i < n_eph; ++i){
    lpd_mk[i] = calc_lpd(oph_mk[i], eph_mk[i], at, beta);
  }
  return(lpd_mk);
}

// [[Rcpp::export]]
std::vector<double> calc_lpd_gtc(std::vector<double> oph_mk,
                                 std::vector<double> eph_mk,
                                 int at,
                                 double beta,
                                 int n_peak,
                                 int n_gtc,
                                 std::vector<int> index_dedup_2_mk){
  std::vector<double> lpd_mk = calc_lpd_mk(oph_mk, eph_mk, at, beta);
  std::vector<double> lpd_gtc(n_gtc, 0);

  int index_start = 0;
  for(int i = 0; i < n_gtc; ++i){
    std::vector<int> index_dedup_2_gtc = extract_subvec_int(index_dedup_2_mk, index_start, index_start + n_peak - 1);
    for(int j = 0; j < n_peak; ++j){
      lpd_gtc[i] += lpd_mk[index_dedup_2_gtc[j]];
    }
    index_start += n_peak;
  }
  return(lpd_gtc);
}

// [[Rcpp::export]]
std::vector<std::vector<double>> calc_lpd_gtc_all(int n_mk,
                                                  std::vector<int> n_dedup_1_mk,
                                                  std::vector<double> oph_dedup_1,
                                                  std::vector<double> eph_dedup_1,
                                                  std::vector<int> at_mk,
                                                  double beta,
                                                  std::vector<int> n_peak_mk,
                                                  std::vector<int> n_gtc_mk,
                                                  std::vector<std::vector<int>> index_dedup_2){
  std::vector<std::vector<double>> lpd_gtc_all(n_mk, std::vector<double>());
  int eph_start = 0;
  for(int i = 0; i < n_mk; ++i){
    int n_dedup_eph = n_dedup_1_mk[i];
    std::vector<double> oph_mk = extract_subvec_double(oph_dedup_1, eph_start, eph_start + n_dedup_eph - 1);
    std::vector<double> eph_mk = extract_subvec_double(eph_dedup_1, eph_start, eph_start + n_dedup_eph - 1);
    int at = at_mk[i];
    int n_peak = n_peak_mk[i];
    int n_gtc = n_gtc_mk[i];

    std::vector<int> index_dedup_2_mk = index_dedup_2.at(i);
    lpd_gtc_all.at(i) = calc_lpd_gtc(oph_mk, eph_mk, at, beta, n_peak, n_gtc, index_dedup_2_mk);

    eph_start += n_dedup_eph;
  }
  return(lpd_gtc_all);
}
