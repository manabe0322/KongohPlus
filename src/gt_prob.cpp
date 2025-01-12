#include "header.h"

// [[Rcpp::export]]
double calc_gtc_prob(std::vector<double> gtc,
                     std::vector<double> pop_al,
                     std::vector<double> pop_prob,
                     double theta,
                     std::vector<double> prof_ref_99,
                     std::vector<int> hyp_id,
                     std::vector<int> n_obs,
                     std::vector<double> obs_al){
  int hnc = hyp_id.size();
  double gtc_prob = 1;

  for(int i = 0; i < hnc; ++i){
    int id = hyp_id[i];
    double al_1 = gtc[2 * i];
    double al_2 = gtc[2 * i + 1];
    if(id != 0){
      if(al_1 != prof_ref_99[2 * id - 2] || al_2 != prof_ref_99[2 * id - 1]){
        gtc_prob *= 0;
        break;
      }
    }else{
      std::vector<double> gt(2);
      gt[0] = al_1;
      gt[1] = al_2;

      double gt_prob = 1;
      int sum_n_obs = accumulate(n_obs.begin(), n_obs.end(), 0);

      for(int j = 0; j < 2; ++j){
        double al = gt[j];
        int index_pop_al = find_index_double(al, pop_al);
        double al_prob = pop_prob[index_pop_al];
        int index_obs_al = find_index_double(al, obs_al);
        int n_al = n_obs[index_obs_al];
        gt_prob *= (n_al * theta + al_prob * (1 - theta)) / (1 + (sum_n_obs - 1) * theta);
        n_obs[index_obs_al] += 1;
        sum_n_obs += 1;
      }

      if(gt[0] != gt[1]){
        gt_prob *= 2;
      }

      gtc_prob *= gt_prob;
    }
  }

  return(gtc_prob);
}


// [[Rcpp::export]]
std::vector<double> calc_gtc_prob_mk(std::vector<std::vector<double>> gtc_mk,
                                     std::vector<double> pop_al,
                                     std::vector<double> pop_prob,
                                     double theta,
                                     std::vector<double> prof_ref_99,
                                     std::vector<int> hyp_id,
                                     std::vector<int> n_obs,
                                     std::vector<double> obs_al){
  int n_gtc_mk = gtc_mk.size();
  std::vector<double> gtc_prob_mk(n_gtc_mk);

  for(int i = 0; i < n_gtc_mk; ++i){
    std::vector<double> gtc = gtc_mk.at(i);
    gtc_prob_mk[i] = calc_gtc_prob(gtc, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs, obs_al);
  }

  return(gtc_prob_mk);
}
