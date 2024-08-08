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

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV 1.27

const double sigma_0 = 2.99416e+01; //mb
const double Q_0 = 1; //GeV
const double x_0 = 7.67079e-05;
const double lambda_star = 3.64361e-01;

const double r_limit = 34.64; // 34.64
const double b_min_limit = 17.32; // 17.32

const int warmup_calls = 10000;
const int integration_calls = 100000;
const int integration_iterations = 1;

const double max_theta_root_excess = 1e-6;

static const double M_X = 3; 

double dipole_amplitude(double r, double x, double b) {
  return exp(-(2*M_PI*b*b)/sigma_0)*(1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow(x/x_0, lambda_star/2)))));
}

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
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
  return_values[0] = acos((A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[0])) {
    return_values[0] = 10; // 10 means that the theta_bar value is invalid and should be skipped
  }
  return_values[1] = acos((-A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[1])) {
    return_values[1] = 10;
  }
  return_values[2] = -acos((A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[2])) {
    return_values[2] = 10;
  }
  return_values[3] = -acos((-A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[3])) {
    return_values[3] = 10;
  }
  return 0;
}

double calc_b_bar(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar) {
  return 1/cos(theta_bar)*(b_min + r/2*cos(phi) - r_bar/2*cos(theta_bar+phi_bar));
}

double trans_T_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {
  if (z*(1-z)*M_X*M_X-m_f*m_f < 0) {
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
    double sub_integrand = 16/gsl_pow_2(2*M_PI)*alpha_em*N_c*e_f*e_f*r*b_min*r_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*M_X*M_X-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))*z*(1-z)
    *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)+epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*cos(-theta_bar[i]+phi-phi_bar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*r_bar))
    *dipole_amplitude(r, x, sqrt(b_min*b_min+2*(1-z)*b_min*r*cos(phi)+gsl_pow_2(1-z)*r*r))*dipole_amplitude(r_bar, x, sqrt(b_min_bar*b_min_bar+gsl_pow_2(1-z)*r_bar*r_bar+(1-z)*2*r_bar*b_min_bar*cos(phi_bar)));

    if (gsl_isnan(sub_integrand)) {
      cout << "T sub_integrand " << i << " is nan" << endl;

      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  return total_integrand;
}

struct parameters {double Q2; double x; double beta;};

double trans_T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return trans_T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->Q2, par->x, par->beta);
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

void trans_integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 6;
  double res, err;

  double xl[dim] = {0, 0, 0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, r_limit, M_PI, 1};

  struct parameters params = {1, 1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  params.beta = par.beta;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;
  double &sigma_fit = par.sigma_fit;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&trans_T_g, dim, &params};

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

double orig_T_integrand(double rx, double ry, double rx_bar, double ry_bar, double bx, double by, double z, double Q2, double x) {
    if (z*(1-z)*M_X*M_X-m_f*m_f<0) {
        return 0;
    }
    return 1/gsl_pow_3(2*M_PI)*alpha_em*N_c*e_f*e_f
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*M_X*M_X-m_f*m_f)*sqrt(gsl_pow_2(rx-rx_bar)+gsl_pow_2(ry-ry_bar)))*z*(1-z)
    *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*sqrt(rx*rx+ry*ry))*gsl_sf_bessel_K0(epsilon(z, Q2)*sqrt(rx_bar*rx_bar+ry_bar*ry_bar))+epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*(rx*rx_bar+ry*ry_bar)/(sqrt(rx*rx+ry*ry)*sqrt(rx_bar*rx_bar+ry_bar*ry_bar))*gsl_sf_bessel_K1(epsilon(z, Q2)*sqrt(rx*rx+ry*ry))*gsl_sf_bessel_K1(epsilon(z, Q2)*sqrt(rx_bar*rx_bar+ry_bar*ry_bar)))*dipole_amplitude(sqrt(rx*rx+ry*ry), x, sqrt(bx*bx+by*by))*dipole_amplitude(sqrt(rx_bar*rx_bar+ry_bar*ry_bar), x, sqrt(bx*bx+by*by));
}

double orig_T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return orig_T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], k[6], par->Q2, par->x);
}

void orig_integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 7;
  double res, err;

  double xl[dim] = {0, 0, 0, 0, 0, 0, 0};
  double xu[dim] = {r_limit, r_limit, r_limit, r_limit, b_min_limit, b_min_limit, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;
  double &sigma_fit = par.sigma_fit;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&orig_T_g, dim, &params};

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

  double Q2 = 1;
  double x = 1e-4;
  double beta = -1;
  double trans_sigma;
  double trans_sigma_error;
  double trans_sigma_fit;

  thread_par_struct trans_parameters(Q2, x, beta, trans_sigma, trans_sigma_error, trans_sigma_fit);
  trans_integrate_for_T_sigma(trans_parameters);

  cout << "trans T sigma=" << trans_sigma << endl;
  cout << "trans T sigma_error=" << trans_sigma_error << endl;
  cout << "trans T sigma_fit=" << trans_sigma_fit << endl;

  double orig_sigma;
  double orig_sigma_error;
  double orig_sigma_fit;

  thread_par_struct trans_parameters(Q2, x, beta, orig_sigma, orig_sigma_error, orig_sigma_fit);
  trans_integrate_for_T_sigma(trans_parameters);

  cout << "orig T sigma=" << orig_sigma << endl;
  cout << "orig T sigma_error=" << orig_sigma_error << endl;
  cout << "orig T sigma_fit=" << orig_sigma_fit << endl;
}