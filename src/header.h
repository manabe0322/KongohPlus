#include <Rcpp.h>
#include <iostream>   // push_back
#include <string>
#include <vector>
#include <boost/math/distributions/gamma.hpp>
using namespace Rcpp;

std::vector<int> extract_subvec_int(std::vector<int> vec, int index_start, int index_end);
std::vector<double> extract_subvec_double(std::vector<double> vec, int index_start, int index_end);
int find_index_double(double element, std::vector<double> vec);
bool is_element_all_vec_double(std::vector<double> vec_a, std::vector<double> vec_b);
std::vector<double> setdiff_double(std::vector<double> vec_num1, std::vector<double> vec_num2);
std::vector<int> setdiff_int(std::vector<int> vec_int1, std::vector<int> vec_int2);

std::vector<double> calc_bs(std::vector<double> alleles);
std::vector<double> calc_fs(std::vector<double> alleles);
std::vector<double> calc_ds(std::vector<double> alleles);
std::vector<double> calc_m2s(std::vector<double> alleles);
std::vector<double> make_cand_peak(std::vector<double> alleles, std::vector<bool> sr_consider);
