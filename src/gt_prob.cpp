#include "header.h"

// [[Rcpp::export]]
std::vector<int> make_n_obs(std::vector<double> pop_al, std::vector<double> prof_ref){
  int n_ref_al = prof_ref.size();
  int index_99 = pop_al.size() - 1;
  std::vector<int> n_obs(index_99 + 1, 0);

  if(n_ref_al != 0){
    for(int i = 0; i < n_ref_al; ++i){
      int index_al = find_index_double(prof_ref[i], pop_al);
      if(index_al != -1){
        n_obs[index_al] += 1;
      }else{
        n_obs[index_99] += 1;
      }
    }
  }

  return(n_obs);
}

// [[Rcpp::export]]
double calc_gt_prob(std::vector<double> gt,
                    std::vector<double> pop_al,
                    std::vector<double> pop_prob,
                    double theta,
                    std::vector<int> n_obs){
  double gt_prob = 1;
  int sum_n_obs = accumulate(n_obs.begin(), n_obs.end(), 0);

  for(int i = 0; i < 2; ++i){
    double al = gt[i];
    int index_al = find_index_double(al, pop_al);
    double al_prob = pop_prob[index_al];
    int n_al = n_obs[index_al];
    gt_prob *= (n_al * theta + al_prob * (1 - theta)) / (1 + (sum_n_obs - 1) * theta);
    n_obs[index_al] += 1;
    sum_n_obs += 1;
  }

  if(gt[0] != gt[1]){
    gt_prob *= 2;
  }

  return(gt_prob);
}

// [[Rcpp::export]]
double calc_gtc_prob(std::vector<double> gtc,
                     std::vector<double> pop_al,
                     std::vector<double> pop_prob,
                     double theta,
                     std::vector<double> prof_ref,
                     std::vector<int> hyp_id,
                     std::vector<int> n_obs){
  int hnc = hyp_id.size();
  double gtc_prob = 1;

  for(int i = 0; i < hnc; ++i){
    int id = hyp_id[i];
    double al_1 = gtc[2 * i];
    double al_2 = gtc[2 * i + 1];
    if(id != 0){
      if(al_1 != prof_ref[2 * id - 2] || al_2 != prof_ref[2 * id - 1]){
        gtc_prob *= 0;
        break;
      }
    }else{
      std::vector<double> gt(2);
      gt[0] = al_1;
      gt[1] = al_2;
      gtc_prob *= calc_gt_prob(gt, pop_al, pop_prob, theta, n_obs);
    }
  }

  return(gtc_prob);
}

// [[Rcpp::export]]
std::vector<double> calc_gtc_adopt_prob(std::vector<std::vector<double>> gtc_adopt,
                                        std::vector<double> pop_al,
                                        std::vector<double> pop_prob,
                                        double theta,
                                        std::vector<double> prof_ref,
                                        std::vector<int> hyp_id,
                                        std::vector<int> n_obs){
  int n_gtc_adopt = gtc_adopt.size();
  std::vector<double> gtc_adopt_prob(n_gtc_adopt);

  for(int i = 0; i < n_gtc_adopt; ++i){
    std::vector<double> gtc = gtc_adopt.at(i);
    gtc_adopt_prob[i] = calc_gtc_prob(gtc, pop_al, pop_prob, theta, prof_ref, hyp_id, n_obs);
  }

  return(gtc_adopt_prob);
}
