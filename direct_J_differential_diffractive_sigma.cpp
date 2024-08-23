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

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27;

const double normalization = 1.0/gsl_pow_3(2*M_PI)*alpha_em*N_c*e_f*e_f;

static double r_limit; // 34.64
static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

static int warmup_calls;
static int integration_calls;
const int integration_iterations = 1;

//const string filename_end = "_20mil_85-225";//

const int debug_precision = 10;

const string dipole_amp_type = "bk";
const string nucleus_type = "p";
const string filename_end = "";

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double x) {
  if (nucleus_type == "p") {
    return get_p_dipole_amplitude(p_table, r, b_min, phi, x, false);
  } else if (nucleus_type == "Pb") {
    return get_Pb_dipole_amplitude(Pb_table, r, b_min, phi, x, false);
  } else {
    throw 1;
  }
}

double calc_phi(double r1, double r2, double b1, double b2, double z) {
  double phi;
  double raw_phi = atan2(-r2, -r1) - atan2((b2+(1-z)*r2), (b1+(1-z)*r1));
  if (M_PI < raw_phi) {
    phi = 2*M_PI - raw_phi;
  } else if (-M_PI <= raw_phi && raw_phi < 0) {
    phi = -raw_phi;
  } else if (raw_phi < -M_PI) {
    phi = 2*M_PI + raw_phi;
  } else {
    phi = raw_phi;
  }
  return phi;
}

double max_sus_integrand = 0;

double T_integrand(double r1, double r2, double b1, double b2, double r1bar, double r2bar, double z, double Q2, double x, double beta) {
  if (z*(1-z)*Q2*(1/beta-1)-m_f*m_f < 0) {
    cout << "setting zero, z=" << z << endl;
    return 0;
  }
  
  double r = sqrt(r1*r1 + r2*r2);
  
  double rbar = sqrt(r1bar*r1bar + r2bar*r2bar);
  
  double b = sqrt(b1*b1 + b2*b2);

  double x2 = gsl_pow_2(1-z)*r*r + 2*(1-z)*(r1*b1+r2*b2) + b*b;
  double y2 = z*z*r*r + 2*z*(r1*b1+r2*b2) + b*b;

  double x2bar = gsl_pow_2(1-z)*rbar*rbar + 2*(1-z)*(r1bar*b1+r2bar*b2) + b*b;
  double y2bar = z*z*rbar*rbar + 2*z*(r1bar*b1+r2bar*b2) + b*b;

  double bmin, bminbar, phi, phibar;
  
  if (x2 < y2) {
    bmin = sqrt(gsl_pow_2(b1+(1-z)*r1) + gsl_pow_2(b2+(1-z)*r2));
    phi = calc_phi(r1, r2, b1, b2, z);
  } else {
    bmin = sqrt(gsl_pow_2(b1-(1-z)*r1) + gsl_pow_2(b2-(1-z)*r2));
    phi = calc_phi(-r1, -r2, b1, b2, z);
  }

  if (x2bar < y2bar) {
    bminbar = sqrt(gsl_pow_2(b1+(1-z)*r1bar) + gsl_pow_2(b2+(1-z)*r2bar));
    phibar = calc_phi(r1bar, r2bar, b1, b2, z);
  } else {
    bminbar = sqrt(gsl_pow_2(b1-(1-z)*r1bar) + gsl_pow_2(b2-(1-z)*r2bar));
    phibar = calc_phi(-r1bar, -r2bar, b1, b2, z);
  }




  /*
  if (r > r_limit) {
    cout << "r overlflow, r=" << r << endl;
  }

  if (bmin > b_min_limit) {
    cout << "bmin overlflow, bmin=" << bmin << endl;
  }
  */

  double integrand =  gsl_sf_bessel_J0(sqrt(z*(1-z)*Q2*(1/beta-1)-m_f*m_f)*sqrt(gsl_pow_2(r1-r1bar)+gsl_pow_2(r2-r2bar)))
  *z*(1-z)
  *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*rbar) + epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*(r1*r1bar+r2*r2bar)/(r*rbar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*rbar))
  *dipole_amplitude(r, bmin, phi, x/beta)*dipole_amplitude(rbar, bminbar, phibar, x/beta);

  if (r > r_limit||rbar > r_limit||bmin > b_min_limit||bminbar > b_min_limit) {
    if (integrand > max_sus_integrand) {
      max_sus_integrand = integrand;
    }
  }
  return integrand;
}

struct parameters {double Q2; double x; double beta;};

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], k[6], par->Q2, par->x, par->beta);
}

struct thread_par_struct
{
  double Q2;
  double x;
  double beta;
  double &sigma;
  double &sigma_error;
  double &sigma_fit;
  int seed;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5, double &a6, int a7) : Q2(a1), x(a2), beta(a3), sigma(a4), sigma_error(a5), sigma_fit(a6), seed(a7) {}
};

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 7;
  double res, err;

  //double L_integrand(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {

  double z_min = (1-sqrt(1-4*m_f*m_f/(par.Q2*(1/par.beta-1))))/2;
  double z_max = (1+sqrt(1-4*m_f*m_f/(par.Q2*(1/par.beta-1))))/2;

  double xl[dim] = {-r_limit, -r_limit, -b_min_limit, -b_min_limit, -r_limit, -r_limit, z_min};
  double xu[dim] = {r_limit, r_limit, b_min_limit, b_min_limit, r_limit, r_limit, z_max};

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
  gsl_rng_set(rng, par.seed);

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
  //cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
  if (nucleus_type == "p") {
    load_p_dipole_amplitudes(p_table, filename);
  } else if (nucleus_type == "Pb") {
    load_Pb_dipole_amplitudes(Pb_table, filename);
  } else {
    throw 1;
  }

  if (nucleus_type == "Pb") {
    r_limit = 657; // 34.64, 657
    b_min_limit = 328; // 17.32, 328
    warmup_calls = 100000;
    integration_calls = 5000000;
  } else if (nucleus_type == "p") {
    r_limit = 34.64;
    b_min_limit = 17.32;
    warmup_calls = 100000;
    integration_calls = 100000;
  } else {
    cout << "invalid nucleus type" << endl;
    throw 1;
  }

  double Q2 = 18;
  double beta = 0.4;
  double x = 1e-2;

  double z_min = (1-sqrt(1-4*m_f*m_f/(Q2*(1/beta-1))))/2;
  double z_max = (1+sqrt(1-4*m_f*m_f/(Q2*(1/beta-1))))/2;

  cout << "Q2=" << Q2 << endl;
  cout << "x=" << x << endl;
  cout << "beta=" << beta << endl;
  cout << z_min << " < z < " << z_max << endl;

  int thread_count = 5;
  double sigma[thread_count];
  double sigma_error[thread_count];
  double sigma_fit[thread_count];
  thread threads[thread_count];

  for (int i=0; i < thread_count; i++){
    thread_par_struct parameters(Q2, x, beta, sigma[i], sigma_error[i], sigma_fit[i], i);
    threads[i] = thread(integrate_for_T_sigma, parameters);
  }
  
  for (int i=0; i<thread_count; i++) {
    threads[i].join();
    cout << "seed=" << i << endl;
    cout << "T sigma=" << sigma[i] << endl;
    cout << "T sigma_error=" << sigma_error[i] << endl;
    cout << "T sigma_fit=" << sigma_fit[i] << endl << endl;
  }
  cout << "max_sus_integrand=" << max_sus_integrand << endl;
}
