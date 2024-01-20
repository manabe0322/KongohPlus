#include <Rcpp.h>
#include <iostream>   // push_back
#include <string>
#include <vector>
using namespace Rcpp;

int find_index_double(double element, std::vector<double> vec);
std::vector<double> calc_bs(std::vector<double> alleles);
std::vector<double> calc_fs(std::vector<double> alleles);
std::vector<double> calc_ds(std::vector<double> alleles);
std::vector<double> calc_m2s(std::vector<double> alleles);
std::vector<double> make_cand_peak(std::vector<double> alleles, std::vector<bool> sr_consider);
