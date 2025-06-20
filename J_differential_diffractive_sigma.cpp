#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>
#include <thread>
#include <chrono>
#include <fstream>
#include <ostream>
#include <sstream>
using namespace std;

#include "direct_dipole_amp_reader.h"

//THIS ONLY WORKS FOR PROTONS

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27;

const double normalization = 16.0/gsl_pow_2(2*M_PI)*alpha_em*N_c*e_f*e_f;

const double r_limit = 34.64; // 34.64
const double b_min_limit = 17.32; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

const int warmup_calls = 100000;
const int integration_calls = 100000000;//20 000 000
const int integration_iterations = 1;

//const string filename_end = "_20mil_85-225";//

const int debug_precision = 10;
const double max_theta_root_excess = 1e-6;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

void read_data_file(string filename, vector<double> &Q2_values, vector<double> &beta_values, vector<double> &x_values, vector<double> &x_pom_F2_values, vector<double> &delta_stat_values, vector<double> &delta_sys_values) {
  ifstream data_file(filename);

  cout << "Reading: " << filename << endl;
  string line;
  while(getline (data_file, line)) {

    long unsigned int i = 0;
    string value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    Q2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    beta_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    x_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    x_pom_F2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    delta_stat_values.push_back(stod(value));
    i++;

    value = "";
    while(i < line.size()) {
      value += line[i];
      i++;
    }
    delta_sys_values.push_back(stod(value));
    i++;

  }
  cout << "Finished reading file" << endl;
}

double calc_h(double r, double b_min, double phi) {
  return 4*b_min*b_min + 4*b_min*r*cos(phi) + r*r;
}

double calc_b1(double r, double b_min, double phi) {
  return b_min + r/2*cos(phi);
}

double calc_b2(double r, double b_min, double phi) {
  return r/2*sin(phi);
}

double calc_j(double b2, double r_bar, double phi_bar) {
  return 4*b2*r_bar*sin(phi_bar);
}

double calc_A(double j, double h, double b1, double b2) {
  return sqrt(j*j + h*(16*b1*b1-gsl_pow_2(j/(2*b2))));
}

bool theta_root_invalid(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar) {
  double excess = abs(tan(theta_bar)*(b_min+r/2*cos(phi)-r_bar/2*cos(phi_bar+theta_bar)) + r_bar/2*sin(phi_bar+theta_bar) - r/2*sin(phi));
  if (excess > max_theta_root_excess) {
    return 1;
  } else {
    return 0;
  }
}

int calc_theta_bar(double return_values[4], double r, double b_min, double phi, double r_bar, double phi_bar) {
  double b1 = calc_b1(r, b_min, phi);
  double b2 = calc_b2(r, b_min, phi);
  double h = calc_h(r, b_min, phi);
  if (r_bar*r_bar > (4*h*b1*b1)/(gsl_pow_2(sin(phi_bar))*(h-4*b2*b2))) {
    //r_bar is too large
    return 1;
  }
  double j = calc_j(b2, r_bar, phi_bar);
  double A = calc_A(j, h, b1, b2);
  if (gsl_isnan(A)) {
    cout << "A is nan!" << endl;
  }
  double acos_arg_plus_type = (A+j)/(2*h);
  if (acos_arg_plus_type > 1) {
    cout << "acos arg of type plus changed: " << acos_arg_plus_type << "->1" << endl;
    acos_arg_plus_type = 1;
  } else if (acos_arg_plus_type < -1) {
    cout << "acos arg of type plus changed: " << acos_arg_plus_type << "->-1" << endl;
    acos_arg_plus_type = -1;
  }
  double acos_arg_minus_type = (-A+j)/(2*h);
  if (acos_arg_minus_type > 1) {
    cout << "acos arg of type minus changed: " << acos_arg_minus_type << "->1" << endl;
    acos_arg_minus_type = 1;
  } else if (acos_arg_minus_type < -1) {
    cout << "acos arg of type minus changed: " << acos_arg_minus_type << "->-1" << endl;
    acos_arg_minus_type = -1;
  }
  return_values[0] = acos(acos_arg_plus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[0])) {
    return_values[0] = 10; // 10 means that the theta_bar value is invalid and should be skipped
  }
  return_values[1] = acos(acos_arg_minus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[1])) {
    return_values[1] = 10;
  }
  return_values[2] = -acos(acos_arg_plus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[2])) {
    return_values[2] = 10;
  }
  return_values[3] = -acos(acos_arg_minus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[3])) {
    return_values[3] = 10;
  }
  return 0;
}

double calc_b_bar(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar) {
  return 1/cos(theta_bar)*(b_min + r/2*cos(phi) - r_bar/2*cos(theta_bar+phi_bar));
}

double step(double argument) {
  if (argument >= 0) {
    return 1;
  } else {
    return 0;
  }
}

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double x) {
  return get_p_dipole_amplitude(table, r, b_min, phi, x);
}

double L_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {
  //static auto t1 = chrono::high_resolution_clock::now();
  if (z*(1-z)*Q2*(1/beta-1)-m_f*m_f < 0) {
    return 0;
  }
  double total_integrand = 0;
  double theta_bar[4];
  int r_bar_too_large = calc_theta_bar(theta_bar, r, b_min, phi, r_bar, phi_bar);
  if (r_bar_too_large == 1) {
    return 0;
  }
  for (int i=0; i<4; i++) {
    if (theta_bar[i] == 10) {
      continue;
    }
    double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i]);
    double sub_integrand = r*b_min*r_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*Q2*(1/beta-1)-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))
    *z*(1-z)*4*Q2*z*z*gsl_pow_2(1-z)*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)
    *dipole_amplitude(r, b_min, phi, x)*dipole_amplitude(r_bar, b_min_bar, phi_bar, x);
    if (gsl_isnan(sub_integrand)) {
      cout << "L sub_integrand " << i << " is nan" << endl;
      /*
      cout << "r=" << setprecision(debug_precision) << r << endl;
      cout << "b_min=" << setprecision(debug_precision) << b_min << endl;
      cout << "phi=" << setprecision(debug_precision) << phi << endl;
      cout << "r_bar=" << setprecision(debug_precision) << r_bar << endl;
      cout << "phi_bar=" << setprecision(debug_precision) << phi_bar << endl;
      cout << "z=" << setprecision(debug_precision) << z << endl;
      cout << "Q2=" << setprecision(debug_precision) << Q2 << endl;
      cout << "x=" << setprecision(debug_precision) << x << endl;
      cout << "beta=" << setprecision(debug_precision) << beta << endl;
      cout << "sub_integrand=" << setprecision(debug_precision) << sub_integrand << endl;
      */
      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  //static auto t2 = chrono::high_resolution_clock::now();
  //auto duration = chrono::duration_cast<chrono::nanoseconds>(t2-t1);
  //cout << duration.count() << endl;

  return total_integrand;
}

double T_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {
  //static auto t1 = chrono::high_resolution_clock::now();
  if (z*(1-z)*Q2*(1/beta-1)-m_f*m_f < 0) {
    return 0;
  }
  double total_integrand = 0;
  double theta_bar[4];
  int r_bar_too_large = calc_theta_bar(theta_bar, r, b_min, phi, r_bar, phi_bar);
  if (r_bar_too_large == 1) {
    return 0;
  }
  for (int i=0; i<4; i++) {
    if (theta_bar[i] == 10) {
      continue;
    }
    double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i]);
    double sub_integrand = r*b_min*r_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*Q2*(1/beta-1)-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))*z*(1-z)
    *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)+epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*cos(-theta_bar[i]+phi-phi_bar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*r_bar))
    *get_p_dipole_amplitude(table, r, b_min, phi, x)*get_p_dipole_amplitude(table, r_bar, b_min_bar, phi_bar, x);

    if (gsl_isnan(sub_integrand)) {
      cout << "T sub_integrand " << i << " is nan" << endl;
      /*
      cout << "r=" << setprecision(debug_precision) << r << endl;
      cout << "b_min=" << setprecision(debug_precision) << b_min << endl;
      cout << "phi=" << setprecision(debug_precision) << phi << endl;
      cout << "r_bar=" << setprecision(debug_precision) << r_bar << endl;
      cout << "phi_bar=" << setprecision(debug_precision) << phi_bar << endl;
      cout << "z=" << setprecision(debug_precision) << z << endl;
      cout << "Q2=" << setprecision(debug_precision) << Q2 << endl;
      cout << "x=" << setprecision(debug_precision) << x << endl;
      cout << "beta=" << setprecision(debug_precision) << beta << endl;
      cout << "sub_integrand=" << setprecision(debug_precision) << sub_integrand << endl;
      */
      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  //static auto t2 = chrono::high_resolution_clock::now();
  //auto duration = chrono::duration_cast<chrono::nanoseconds>(t2-t1);
  //cout << duration.count() << endl;

  return total_integrand;
}

struct parameters {double Q2; double x; double beta;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->Q2, par->x, par->beta);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->Q2, par->x, par->beta);
}

struct thread_par_struct
{
  double Q2;
  double x;
  double beta;
  double &sigma;
  double &sigma_error;
  double &sigma_fit;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5, double &a6) : Q2(a1), x(a2), beta(a3), sigma(a4), sigma_error(a5), sigma_fit(a6) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  const int dim = 6;
  double res, err;

  //double L_integrand(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {
  double xl[6] = {0, 0, 0, 0, 0, 0};
  double xu[6] = {r_limit, b_min_limit, M_PI, r_limit, M_PI, 1};

  struct parameters params = {1, 1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  params.beta = par.beta;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;
  double &sigma_fit = par.sigma_fit;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
  if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    static auto t1 = chrono::high_resolution_clock::now();
    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
    if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    //cout << "L iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(L_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  sigma_fit = gsl_monte_vegas_chisq(L_s);
  cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 6;
  double res, err;

  //double L_integrand(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {
  double xl[6] = {0, 0, 0, 0, 0, 0};
  double xu[6] = {r_limit, b_min_limit, M_PI, r_limit, M_PI, 1};

  struct parameters params = {1, 1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  params.beta = par.beta;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;
  double &sigma_fit = par.sigma_fit;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    static auto t1 = chrono::high_resolution_clock::now();
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    //cout << "T iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(T_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  sigma_fit = gsl_monte_vegas_chisq(T_s);
  cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  string filename = "data/dipole_amplitude_with_IP_dependence_bk_p_diffraction.csv";
  load_p_dipole_amplitudes(table, filename);

  double Q2 = 35;
  double beta = 0.04;
  double x = 0.018*beta;
  double sigma;
  double sigma_error;
  double sigma_fit;

  cout << "Q2=" << Q2 << endl;
  cout << "x=" << x << endl;
  cout << "beta=" << beta << endl;

  thread_par_struct parameters(Q2, x, beta, sigma, sigma_error, sigma_fit);
  integrate_for_T_sigma(parameters);

  cout << "Q2=" << Q2 << endl;
  cout << "x=" << x << endl;
  cout << "beta=" << beta << endl;

  cout << endl;

  cout << "T sigma=" << sigma << endl;
  cout << "T sigma_error=" << sigma_error << endl;
  cout << "T sigma_fit=" << sigma_fit << endl;
}
