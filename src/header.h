#include <Rcpp.h>
#include <iostream>   // push_back
#include <math.h>     // fmod
#include <string>
#include <vector>
using namespace Rcpp;

std::vector<double> calc_bs(std::vector<double> alleles);
std::vector<double> calc_fs(std::vector<double> alleles);
std::vector<double> calc_ds(std::vector<double> alleles);
std::vector<double> calc_m2s(std::vector<double> alleles);
std::vector<double> make_cand_peak(std::vector<double> alleles,
                                   std::vector<bool> sr_consider);
