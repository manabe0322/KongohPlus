#include <Rcpp.h>
#include <math.h>
#include <string>
#include <vector>
using namespace Rcpp;

std::vector<double> calc_bs(std::vector<double> alleles);
std::vector<double> calc_fs(std::vector<double> alleles);
std::vector<double> calc_ds(std::vector<double> alleles);
std::vector<double> calc_m2s(std::vector<double> alleles);
