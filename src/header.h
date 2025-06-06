#include <Rcpp.h>
#include <algorithm>  // max_element
#include <cmath>      // round, pow
#include <iostream>   // push_back
#include <numeric>    // accumulate
#include <omp.h>      // OpenMP
#include <string>
#include <vector>
#include <boost/math/distributions/gamma.hpp>
using namespace Rcpp;

std::vector<int> extract_subvec_int(std::vector<int> vec, int index_start, int index_end);
std::vector<double> extract_subvec_double(std::vector<double> vec, int index_start, int index_end);
int find_index_double(double element, std::vector<double> vec);
std::vector<double> setdiff_double(std::vector<double> vec_num1, std::vector<double> vec_num2);
std::vector<int> setdiff_int(std::vector<int> vec_int1, std::vector<int> vec_int2);

double calc_bs_one(double al);
std::vector<double> calc_bs(std::vector<double> alleles);
double calc_fs_one(double al);
std::vector<double> calc_fs(std::vector<double> alleles);
double calc_ds_one(double al);
std::vector<double> calc_ds(std::vector<double> alleles);
double calc_m2s_one(double al);
std::vector<double> calc_m2s(std::vector<double> alleles);
double calc_parent_al_bs_one(double bs);
double calc_parent_al_fs_one(double fs);
double calc_parent_al_ds_one(double ds);
double calc_parent_al_m2s_one(double m2s);
std::vector<double> extract_gtc(std::vector<std::vector<double>> cand_gts, std::vector<int> ind_gtc, int hnc);
double log_sum_exp(const std::vector<double>& log_values);
