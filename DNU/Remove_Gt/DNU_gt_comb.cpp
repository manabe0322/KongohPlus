bool judge_gtc_phr(std::vector<double> gtc,
                   std::vector<double> peaks,
                   std::vector<double> heights,
                   double phr_fltr,
                   double st){
  bool adopt_gtc = true;
  int hnc = gtc.size() / 2;
  std::vector<int> index_gtc_all(2 * hnc);
  std::iota(index_gtc_all.begin(), index_gtc_all.end(), 0);

  for(int i = 0; i < hnc; ++i){
    double al_1 = gtc[2 * i];
    double al_2 = gtc[2 * i + 1];

    if(al_1 != al_2){
      std::vector<int> index_gtc_current(2);
      index_gtc_current[0] = 2 * i;
      index_gtc_current[1] = 2 * i + 1;
      std::vector<int> index_gtc_other = setdiff_int(index_gtc_all, index_gtc_current);

      int n_gtc_other = 2 * hnc - 2;
      std::vector<double> gtc_other(n_gtc_other);
      for(int j = 0; j < n_gtc_other; ++j){
        gtc_other[j] = gtc[index_gtc_other[j]];
      }

      int index_gtc_other_1 = find_index_double(al_1, gtc_other);
      int index_gtc_other_2 = find_index_double(al_2, gtc_other);
      if(index_gtc_other_1 == -1 && index_gtc_other_2 == -1){
        int index_peak_1 = find_index_double(al_1, peaks);
        int index_peak_2 = find_index_double(al_2, peaks);

        double h1 = heights[index_peak_1];
        double h2 = heights[index_peak_2];
        if(h1 >= st && h2 >= st){
          double phr = h1 / h2;
          if(phr < phr_fltr || phr > 1 / phr_fltr){
            adopt_gtc = false;
            break;
          }
        }
      }
    }
  }
  return(adopt_gtc);
}

// [[Rcpp::export]]
bool judge_gtc_mr(std::vector<double> gtc,
                  std::vector<double> peaks,
                  std::vector<double> heights,
                  double mr_fltr,
                  double st){
  bool adopt_gtc = true;
  int hnc = gtc.size() / 2;
  std::vector<int> index_gtc_all(2 * hnc);
  std::iota(index_gtc_all.begin(), index_gtc_all.end(), 0);
  std::vector<double> uniq_heights;
  std::vector<bool> pass_st;

  for(int i = 0; i < hnc; ++i){
    double al_1 = gtc[2 * i];
    double al_2 = gtc[2 * i + 1];

    std::vector<int> index_gtc_current(2);
    index_gtc_current[0] = 2 * i;
    index_gtc_current[1] = 2 * i + 1;
    std::vector<int> index_gtc_other = setdiff_int(index_gtc_all, index_gtc_current);

    int n_gtc_other = 2 * hnc - 2;
    std::vector<double> gtc_other(n_gtc_other);
    for(int j = 0; j < n_gtc_other; ++j){
      gtc_other[j] = gtc[index_gtc_other[j]];
    }

    int index_gtc_other_1 = find_index_double(al_1, gtc_other);
    int index_gtc_other_2 = find_index_double(al_2, gtc_other);

    double n_divide = 0;
    std::vector<double> uniq_h_per_al(2);
    if(index_gtc_other_1 == -1){
      int index_peak_1 = find_index_double(al_1, peaks);
      uniq_h_per_al[0] = heights[index_peak_1];
      n_divide += 1;
    }
    if(index_gtc_other_2 == -1){
      if(al_1 != al_2){
        int index_peak_2 = find_index_double(al_2, peaks);
        uniq_h_per_al[1] = heights[index_peak_2];
      }
      n_divide += 1;
    }

    if(n_divide != 0){
      double uniq_h = accumulate(uniq_h_per_al.begin(), uniq_h_per_al.end(), 0.0) / n_divide;
      uniq_heights.push_back(uniq_h);
      pass_st.push_back(uniq_h >= st);
    }
  }

  int n_uniq = uniq_heights.size();
  if(n_uniq >= 2){
    for(int i = 0; i < (n_uniq - 1); ++i){
      if(pass_st[i]){
        double uniq_heights_target = uniq_heights[i];
        std::vector<double> uniq_heights_others = extract_subvec_double(uniq_heights, i + 1, n_uniq - 1);
        int n_others = uniq_heights_others.size();
        for(int j = 0; j < n_others; ++j){
          double mr = uniq_heights_target / uniq_heights_others[j];
          if(mr < mr_fltr){
            adopt_gtc = false;
            break;
          }
        }
        if(!adopt_gtc){
          break;
        }
      }
    }
  }
  return(adopt_gtc);
}

// [[Rcpp::export]]
bool is_element_all_vec_double(std::vector<double> vec_a, std::vector<double> vec_b){
  bool ans = true;
  int n_vec_a = vec_a.size();
  for(int i = 0; i < n_vec_a; ++i){
    double element = vec_a[i];
    auto index = std::find(vec_b.begin(), vec_b.end(), element);
    ans = index != vec_b.end();
    if(!ans){
      break;
    }
  }
  return(ans);
}
