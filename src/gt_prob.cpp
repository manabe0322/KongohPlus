#include "header.h"

// [[Rcpp::export]]
double calc_gts_prob(std::vector<double> gts,
                     std::vector<double> pop_al,
                     std::vector<double> pop_prob,
                     double theta,
                     std::vector<double> prof_ref_99,
                     std::vector<int> hyp_id,
                     std::vector<int> n_obs,
                     std::vector<double> obs_al,
                     std::vector<int> index_dropout_cpp){
  int hnc = hyp_id.size();
  double gts_prob = 1;

  for(int i = 0; i < hnc; ++i){
    int id = hyp_id[i];
    double al_1 = gts[2 * i];
    double al_2 = gts[2 * i + 1];
    if(id != 0){
      if(al_1 != prof_ref_99[2 * id - 2] || al_2 != prof_ref_99[2 * id - 1]){
        gts_prob *= 0;
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

        if(!index_dropout_cpp.empty()){
          int n_dropout_al = index_dropout_cpp.size();
          for(int k = 0; k < n_dropout_al; ++k){
            n_al += n_obs[index_dropout_cpp[k]];
          }
        }

        gt_prob *= (n_al * theta + al_prob * (1 - theta)) / (1 + (sum_n_obs - 1) * theta);
        n_obs[index_obs_al] += 1;
        sum_n_obs += 1;
      }

      if(gt[0] != gt[1]){
        gt_prob *= 2;
      }

      gts_prob *= gt_prob;
    }
  }

  return(gts_prob);
}


// [[Rcpp::export]]
std::vector<double> calc_gts_prob_mk(std::vector<std::vector<double>> gts_mk,
                                    std::vector<double> pop_al,
                                    std::vector<double> pop_prob,
                                    double theta,
                                    std::vector<double> prof_ref_99,
                                    std::vector<int> hyp_id,
                                    std::vector<int> n_obs,
                                    std::vector<double> obs_al,
                                    std::vector<int> index_dropout_cpp){
  int n_gts_mk = gts_mk.size();
  std::vector<double> gts_prob_mk(n_gts_mk);

  #pragma omp parallel for
  for(int i = 0; i < n_gts_mk; ++i){
    std::vector<double> gts = gts_mk[i];
    std::vector<int> n_obs_copy = n_obs;
    gts_prob_mk[i] = calc_gts_prob(gts, pop_al, pop_prob, theta, prof_ref_99, hyp_id, n_obs_copy, obs_al, index_dropout_cpp);
  }

  return(gts_prob_mk);
}
