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
#include <iomanip>
using namespace std;

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV 1.27

const double sigma_0 = 2.99416e+01; //mb
const double Q_0 = 1; //GeV
const double x_0 = 7.67079e-05;
const double lambda_star = 3.64361e-01;

const double r_limit = 50; // 34.64
const double b_min_limit = 50; // 17.32

const int warmup_calls = 100000;
const int integration_calls = 10000000;
const int integration_iterations = 1;

const double max_theta_root_excess = 1e-6;
const int debug_precision = 10;
const double value_for_nan = 1e-5;
const double max_value = 1e10;

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

double calc_h(double r, double b_min, double phi, double z) {
  return b_min*b_min + (1-z)*2*b_min*r*cos(phi) + gsl_pow_2(1-z)*r*r;
}

double calc_b1(double r, double b_min, double phi, double z) {
  return b_min + r*(1-z)*cos(phi);
}

double calc_b2(double r, double b_min, double phi, double z) {
  return r*(1-z)*sin(phi);
}

double calc_j(double b2, double r_bar, double phi_bar, double z) {
  return (1-z)*2*b2*r_bar*sin(phi_bar);
}

double calc_A(double j, double h, double b1, double b2) {
  return sqrt(j*j + 4*h*(b1*b1-gsl_pow_2(j/(2*b2))));
}

bool theta_root_invalid(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar, double z) {
  double excess = abs(tan(theta_bar)*(b_min+r*(1-z)*cos(phi)-r_bar*(1-z)*cos(phi_bar+theta_bar)) + r_bar*(1-z)*sin(phi_bar+theta_bar) - r*(1-z)*sin(phi));
  if (excess > max_theta_root_excess) {
    return 1;
  } else {
    return 0;
  }
}

int calc_theta_bar(double return_values[4], double r, double b_min, double phi, double r_bar, double phi_bar, double z) {
  double b1 = calc_b1(r, b_min, phi, z);
  double b2 = calc_b2(r, b_min, phi, z);
  double h = calc_h(r, b_min, phi, z);
  if (r_bar*r_bar > (4*h*b1*b1)/(gsl_pow_2((1-z)*2*sin(phi_bar))*(h-b2*b2))) {
    //r_bar is too large
    return 1;
  }
  double j = calc_j(b2, r_bar, phi_bar, z);
  double A = calc_A(j, h, b1, b2);
  if (gsl_isnan(A)) {
    cout << "A is nan!" << endl;
    cout << "j=" << j << endl;
    cout << "h=" << h << endl;
    cout << "b1=" << b1 << endl;
    cout << "b2=" << b2 << endl;
    cout << "A'=" << j*j + 4*h*(b1*b1-gsl_pow_2(j/(2*b2))) << endl;
    cout << endl;

  }
  double acos_arg_plus_type = (A+j)/(2*h);
  if (acos_arg_plus_type > 1) {
    cout << "acos arg of type plus changed: " << setprecision(debug_precision) << acos_arg_plus_type << "->1" << endl;
    acos_arg_plus_type = 1;
  } else if (acos_arg_plus_type < -1) {
    cout << "acos arg of type plus changed: " << setprecision(debug_precision) << acos_arg_plus_type << "->-1" << endl;
    acos_arg_plus_type = -1;
  }
  double acos_arg_minus_type = (-A+j)/(2*h);
  if (acos_arg_minus_type > 1) {
    cout << "acos arg of type minus changed: " << setprecision(debug_precision) << acos_arg_minus_type << "->1" << endl;
    acos_arg_minus_type = 1;
  } else if (acos_arg_minus_type < -1) {
    cout << "acos arg of type minus changed: " << setprecision(debug_precision) << acos_arg_minus_type << "->-1" << endl;
    acos_arg_minus_type = -1;
  }
  return_values[0] = acos(acos_arg_plus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[0], z)) {
    return_values[0] = 10; // 10 means that the theta_bar value is invalid and should be skipped
  }
  return_values[1] = acos(acos_arg_minus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[1], z)) {
    return_values[1] = 10;
  }
  return_values[2] = -acos(acos_arg_plus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[2], z)) {
    return_values[2] = 10;
  }
  return_values[3] = -acos(acos_arg_minus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[3], z)) {
    return_values[3] = 10;
  }
  return 0;
}

double calc_b_bar(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar, double z) {
  return 1/cos(theta_bar)*(b_min + (1-z)*r*cos(phi) - (1-z)*r_bar*cos(theta_bar+phi_bar));
}

double trans_T_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double Q2, double x, double beta, double theta_selection[4]) {
  double z = 1.0/2;
  if (z*(1-z)*M_X*M_X-m_f*m_f < 0) {
    return 0;
  }
  double total_integrand = 0;
  double theta_bar[4];
  int r_bar_too_large = calc_theta_bar(theta_bar, r, b_min, phi, r_bar, phi_bar, z);
  if (r_bar_too_large == 1) {
    return 0;
  }
  for (int i=0; i<4; i++) {
    if (theta_selection[i] != 1) {
      continue;
    }
    if (theta_bar[i] == 10) {
      continue;
    }
    double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i], z);

    //cout << sqrt(b_min*b_min+2*(1-z)*b_min*r*cos(phi)+gsl_pow_2(1-z)*r*r) - sqrt(b_min_bar*b_min_bar+gsl_pow_2(1-z)*r_bar*r_bar+(1-z)*2*r_bar*b_min_bar*cos(phi_bar)) << endl;
    //cout << "r=" << r << ", bmin=" << b_min << ", phi=" << phi << "r_bar=" << r_bar << ", phi_bar=" << phi_bar << ", z=" << z << ", Q2=" << Q2 << ", x=" << x << ", beta=" << beta << ", M_X=" << M_X << endl;
    //cout << "b=" << sqrt(b_min*b_min+2*(1-z)*b_min*r*cos(phi)+gsl_pow_2(1-z)*r*r) << ", b_bar=" << sqrt(b_min_bar*b_min_bar+gsl_pow_2(1-z)*r_bar*r_bar+(1-z)*2*r_bar*b_min_bar*cos(phi_bar)) << endl << endl;
    
    double sub_integrand = 1/abs(cos(theta_bar[i]))*1/abs(b_min_bar*cos(theta_bar[i])+(1-z)*r_bar*cos(theta_bar[i]+phi_bar))
    *16/gsl_pow_2(2*M_PI)*alpha_em*N_c*e_f*e_f*r*b_min*r_bar*b_min_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*M_X*M_X-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))*z*(1-z)
    *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)+epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*cos(-theta_bar[i]+phi-phi_bar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*r_bar))
    *dipole_amplitude(r, x, sqrt(b_min*b_min+2*(1-z)*b_min*r*cos(phi)+gsl_pow_2(1-z)*r*r))*dipole_amplitude(r_bar, x, sqrt(b_min_bar*b_min_bar+gsl_pow_2(1-z)*r_bar*r_bar+(1-z)*2*r_bar*b_min_bar*cos(phi_bar)));

    if (gsl_isnan(sub_integrand)) {
      cout << "T sub_integrand " << i << " is nan" << endl;
      cout << "r=" << r << ", bmin=" << b_min << ", phi=" << phi << "r_bar=" << r_bar << ", phi_bar=" << phi_bar << ", z=" << z << ", Q2=" << Q2 << ", x=" << x << ", beta=" << beta << ", M_X=" << M_X << endl << endl;;

      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  return total_integrand;
}

struct parameters {double Q2; double x; double beta; double theta_1; double theta_2; double theta_3; double theta_4;};

double trans_T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  double theta_selection[4] = {par->theta_1, par->theta_2, par->theta_3, par->theta_4};
  return trans_T_integrand(k[0], k[1], k[2], k[3], k[4], par->Q2, par->x, par->beta, theta_selection);
}

struct thread_par_struct
{
  double Q2;
  double x;
  double beta;
  double &sigma;
  double &sigma_error;
  double &sigma_fit;
  int theta_1;
  int theta_2;
  int theta_3;
  int theta_4;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5, double &a6, double a7, double a8, double a9, double a10) : Q2(a1), x(a2), beta(a3), sigma(a4), sigma_error(a5), sigma_fit(a6), theta_1(a7), theta_2(a8), theta_3(a9), theta_4(a10) {}
};

void trans_integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 5;
  double res, err;

  double xl[dim] = {0, 0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, r_limit, M_PI};

  struct parameters params = {1, 1, 1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  params.beta = par.beta;
  params.theta_1 = par.theta_1;
  params.theta_2 = par.theta_2;
  params.theta_3 = par.theta_3;
  params.theta_4 = par.theta_4;
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
    res = value_for_nan;
    cout << "nan found at x=" << params.x << endl;
  }
  if (res > max_value) {
    res = value_for_nan;
    cout << "inf found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  sigma_fit = gsl_monte_vegas_chisq(T_s);
  cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

double orig_T_integrand(double rx, double ry, double rx_bar, double ry_bar, double b, double Q2, double x) {
    double z = 1.0/2;
    if (z*(1-z)*M_X*M_X-m_f*m_f<0) {
        return 0;
    }
    double r = sqrt(rx*rx+ry*ry);
    double r_bar = sqrt(rx_bar*rx_bar+ry_bar*ry_bar);
    return 2*M_PI*b*1/gsl_pow_3(2*M_PI)*alpha_em*N_c*e_f*e_f // 6.7065002*
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*M_X*M_X-m_f*m_f)*sqrt(gsl_pow_2(rx-rx_bar)+gsl_pow_2(ry-ry_bar)))*z*(1-z)
    *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)+epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*(rx*rx_bar+ry*ry_bar)/(r*r_bar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*r_bar))*dipole_amplitude(r, x, b)*dipole_amplitude(r_bar, x, b);
}

double orig_T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return orig_T_integrand(k[0], k[1], k[2], k[3], k[4], par->Q2, par->x);
}

void orig_integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 5;
  double res, err;

  double xl[dim] = {-r_limit, -r_limit, -r_limit, -r_limit, 0};
  double xu[dim] = {r_limit, r_limit, r_limit, r_limit, b_min_limit};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  params.beta = par.beta;
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
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, 1000000, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    //cout << "T iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(T_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = value_for_nan;
    cout << "nan found at x=" << params.x << endl;
  }
  if (res > max_value) {
    res = value_for_nan;
    cout << "inf found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  sigma_fit = gsl_monte_vegas_chisq(T_s);
  cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  //r=6.6823, bmin=13.6696, phi=1.86557r_bar=16.187, phi_bar=1.15127, z=0.683371, Q2=0, x=0.00188739, beta=-1, M_X=3.00891
  //b=13.2109, b_bar=11.2299

  /*
  double results[4];

  double r = 6.6823;
  double b_min = 13.6696;
  double phi = 1.86557;
  double r_bar = 16.187;
  double phi_bar = 1.15127;
  double z = 0.683371;

  calc_theta_bar(results, r, b_min, phi, r_bar, phi_bar, z);

  cout << results[0] << endl;
  cout << results[1] << endl;
  cout << results[2] << endl;
  cout << results[3] << endl;
  cout << "b_min_bar=" << calc_b_bar(r, b_min, phi, r_bar, phi_bar, results[0], z) << endl;

  cout << endl;

  cout << "b1=" << calc_b1(r, b_min, phi, z) << endl;
  cout << "b2=" << calc_b2(r, b_min, phi, z) << endl;
  cout << "j=" << calc_j(calc_b2(r, b_min, phi, z), r_bar, phi_bar, z) << endl;

  return 0;
  */

  gsl_set_error_handler_off();

  const double Q2_values[] = {0}; //{0, 1, 2, 3, 4, 5}

  const int x_steps = 5;
  const double x_start = 1e-5;
  const double x_stop = 0.01;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  double beta = -1;

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Variable change comparison;x;cross section (mb)");
  
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps], T_sigma_fits[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, beta, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i], 1, 1, 1, 1);
      T_threads[i] = thread(trans_integrate_for_T_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, T_x_values, T_sigma_values, T_x_errors, T_sigma_errors);
    TString subgraph_name = "trans Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    subgraph->SetLineColor(j+1);
    T_graphs->Add(subgraph, "PC");
  }

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps], T_sigma_fits[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, beta, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i], 1, 0, 0, 0);
      T_threads[i] = thread(trans_integrate_for_T_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, T_x_values, T_sigma_values, T_x_errors, T_sigma_errors);
    TString subgraph_name = "theta 1 Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    subgraph->SetLineColor(2);
    subgraph->SetLineStyle(3);
    T_graphs->Add(subgraph, "PC");
  }

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps], T_sigma_fits[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, beta, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i], 0, 1, 0, 0);
      T_threads[i] = thread(trans_integrate_for_T_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, T_x_values, T_sigma_values, T_x_errors, T_sigma_errors);
    TString subgraph_name = "theta 2 Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    subgraph->SetLineColor(3);
    subgraph->SetLineStyle(4);
    T_graphs->Add(subgraph, "PC");
  }

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps], T_sigma_fits[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, beta, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i], 0, 0, 1, 0);
      T_threads[i] = thread(trans_integrate_for_T_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, T_x_values, T_sigma_values, T_x_errors, T_sigma_errors);
    TString subgraph_name = "theta 3 Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    subgraph->SetLineColor(4);
    subgraph->SetLineStyle(5);
    T_graphs->Add(subgraph, "PC");
  }

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps], T_sigma_fits[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, beta, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i], 0, 0, 0, 1);
      T_threads[i] = thread(trans_integrate_for_T_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, T_x_values, T_sigma_values, T_x_errors, T_sigma_errors);
    TString subgraph_name = "theta 4 Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    subgraph->SetLineColor(5);
    subgraph->SetLineStyle(6);
    T_graphs->Add(subgraph, "PC");
  }

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps], T_sigma_fits[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, beta, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i], 0, 0, 0, 0);
      T_threads[i] = thread(orig_integrate_for_T_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, T_x_values, T_sigma_values, T_x_errors, T_sigma_errors);
    TString subgraph_name = "orig Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    subgraph->SetLineStyle(2);
    subgraph->SetLineColor(j+1);
    T_graphs->Add(subgraph, "PC");
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  T_graphs->Draw("A");

  gPad->SetLogx();
  gPad->SetLogy();

  comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  comparison_canvas->Print("figures/variable_change_comparison.pdf");

  return 0;
}
