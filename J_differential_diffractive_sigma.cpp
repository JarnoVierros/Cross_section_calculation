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
const double m_f = 1.27; //GeV 1.27

const double normalization = 16/gsl_pow_2(2*M_PI)*alpha_em*N_c*e_f*e_f;

const double r_limit = 4; // 34.64
const double b_min_limit = 10; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

const int warmup_calls = 10000;
const int integration_calls = 1000000;//20000000
const int integration_iterations = 10;

const int debug_precision = 20;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

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
  return_values[1] = acos((-A+j)/(2*h));
  return_values[2] = -acos((A+j)/(2*h));
  return_values[3] = -acos((-A+j)/(2*h));
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
  return get_dipole_amplitude(table, r, b_min, phi, x);
}

double L_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double z, double Q2, double x_pom, double beta) {
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
    double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i]);
    double sub_integrand = r*b_min*r_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*Q2*(1/beta-1)-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))
    *z*(1-z)*4*Q2*z*z*gsl_pow_2(1-z)*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)
    *get_dipole_amplitude(table, r, b_min, phi, beta*x_pom)*get_dipole_amplitude(table, r_bar, b_min_bar, phi_bar, beta*x_pom);

    if (gsl_isnan(sub_integrand)) {
      cout << "L sub_integrand " << i << " is nan" << endl;
      cout << "r=" << setprecision(debug_precision) << r << endl;
      cout << "b_min=" << setprecision(debug_precision) << b_min << endl;
      cout << "phi=" << setprecision(debug_precision) << phi << endl;
      cout << "r_bar=" << setprecision(debug_precision) << r_bar << endl;
      cout << "phi_bar=" << setprecision(debug_precision) << phi_bar << endl;
      cout << "z=" << setprecision(debug_precision) << z << endl;
      cout << "Q2=" << setprecision(debug_precision) << Q2 << endl;
      cout << "x_pom=" << setprecision(debug_precision) << x_pom << endl;
      cout << "beta=" << setprecision(debug_precision) << beta << endl;
      cout << "sub_integrand=" << setprecision(debug_precision) << sub_integrand << endl;
      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  //static auto t2 = chrono::high_resolution_clock::now();
  //auto duration = chrono::duration_cast<chrono::nanoseconds>(t2-t1);
  //cout << duration.count() << endl;

  return total_integrand;
}

double T_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double z, double Q2, double x_pom, double beta) {
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
    double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i]);
    double sub_integrand = r*b_min*r_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*Q2*(1/beta-1)-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))*z*(1-z)
    *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)+epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*cos(-theta_bar[i]+phi-phi_bar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*r_bar))
    *get_dipole_amplitude(table, r, b_min, phi, beta*x_pom)*get_dipole_amplitude(table, r_bar, b_min_bar, phi_bar, beta*x_pom);

    if (gsl_isnan(sub_integrand)) {
      cout << "T sub_integrand " << i << " is nan" << endl;
      cout << "r=" << setprecision(debug_precision) << r << endl;
      cout << "b_min=" << setprecision(debug_precision) << b_min << endl;
      cout << "phi=" << setprecision(debug_precision) << phi << endl;
      cout << "r_bar=" << setprecision(debug_precision) << r_bar << endl;
      cout << "phi_bar=" << setprecision(debug_precision) << phi_bar << endl;
      cout << "z=" << setprecision(debug_precision) << z << endl;
      cout << "Q2=" << setprecision(debug_precision) << Q2 << endl;
      cout << "x_pom=" << setprecision(debug_precision) << x_pom << endl;
      cout << "beta=" << setprecision(debug_precision) << beta << endl;
      cout << "sub_integrand=" << setprecision(debug_precision) << sub_integrand << endl;
      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  //static auto t2 = chrono::high_resolution_clock::now();
  //auto duration = chrono::duration_cast<chrono::nanoseconds>(t2-t1);
  //cout << duration.count() << endl;

  return total_integrand;
}

struct parameters {double Q2; double x_pom; double beta;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->Q2, par->x_pom, par->beta);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->Q2, par->x_pom, par->beta);
}

struct thread_par_struct
{
  double Q2;
  double x_pom;
  double beta;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5) : Q2(a1), x_pom(a2), beta(a3), sigma(a4), sigma_error(a5) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  const int dim = 6;
  double res, err;

  //double L_integrand(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double z, double Q2, double x_pom, double beta) {
  double xl[6] = {0, 0, 0, 0, 0, 0};
  double xu[6] = {r_limit, b_min_limit, M_PI, r_limit, M_PI, 1};

  struct parameters params = {1, 1, 1};
  params.Q2 = par.Q2;
  params.x_pom = par.x_pom;
  params.beta = par.beta;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
  if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    static auto t1 = chrono::high_resolution_clock::now();
    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
    if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    cout << "L iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(L_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x_pom=" << params.x_pom << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "L, Q²=" << params.Q2 << ", x_pom=" << params.x_pom << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 6;
  double res, err;

  //double L_integrand(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double z, double Q2, double x_pom, double beta) {
  double xl[6] = {0, 0, 0, 0, 0, 0};
  double xu[6] = {r_limit, b_min_limit, M_PI, r_limit, M_PI, 1};

  struct parameters params = {1, 1, 1};
  params.Q2 = par.Q2;
  params.x_pom = par.x_pom;
  params.beta = par.beta;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    static auto t1 = chrono::high_resolution_clock::now();
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    cout << "T iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(T_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x_pom=" << params.x_pom << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "T, Q²=" << params.Q2 << ", x_pom=" << params.x_pom << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

/*
  static auto t1 = chrono::high_resolution_clock::now();
  int status = gsl_integration_qags (&F_phi, 0, M_PI, absolute_precision, 0.01, integration_limit, w_phi, &result, &error);
  static auto t2 = chrono::high_resolution_clock::now();
  handle_status(status);
  auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);

  cout << "bmin= " << bmin << ", time=" << duration.count() << endl;
  cout << "average phibar time: " << phibar_t_total/phibar_t_count << ", average rbar time: " << rbar_t_total/rbar_t_count << endl;
  phi_t_total += duration.count();
  phi_t_count++;
*/

/*
void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 4;
  double res, err;

  double xl[4] = {0, 0, 0, 0};
  double xu[4] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}
*/

int main() {

  gsl_set_error_handler_off();

  string filename = "data/dipole_amplitude_with_IP_dependence.csv";
  load_dipole_amplitudes(table, filename);

  double L_sigma;
  double L_sigma_error;
  thread_par_struct L_par(4.5, 0.00012, 0.04, L_sigma, L_sigma_error);
  thread L_integration = thread(integrate_for_L_sigma, L_par);

  double T_sigma;
  double T_sigma_error;
  thread_par_struct T_par(4.5, 0.00012, 0.04, T_sigma, T_sigma_error);
  thread T_integration = thread(integrate_for_T_sigma, T_par);

  L_integration.join();
  T_integration.join();

  cout << "L sigma: " << L_sigma << ", error: " << L_sigma_error << endl;
  cout << "T sigma: " << T_sigma << ", error: " << T_sigma_error << endl;

  /*

  const int Q2_values[] = {1, 3, 5, 8, 10};

  const int x_steps = 30;
  const double x_start = 1e-5;
  const double x_stop = 0.01;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  string filename = "data/dipole_amplitude_with_IP_dependence.csv";
  load_dipole_amplitudes(table, filename);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Diffractive longitudinal cross section;x;cross section (mb)");

  ofstream L_output_file("data/diff_L_sigma_x.txt");
  L_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting L integration" << endl;
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double L_x_values[x_steps], L_sigma_values[x_steps], L_x_errors[x_steps], L_sigma_errors[x_steps];
    thread L_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      L_x_values[i] = x;
      L_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, L_sigma_values[i], L_sigma_errors[i]);
      L_threads[i] = thread(integrate_for_L_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      L_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, L_x_values, L_sigma_values, L_x_errors, L_sigma_errors);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);

    for (int i=0; i<x_steps; i++) {
      ostringstream x;
      x << L_x_values[i];
      ostringstream sigma;
      sigma << L_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << L_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + x.str() + ";" + sigma.str() + ";" + sigma_err.str();
      L_output_file << line << endl;
    }
  }
  L_output_file.close();

  TCanvas* L_sigma_canvas = new TCanvas("diff_L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  L_sigma_canvas->Print("figures/diff_L_sigma_x_distribution.pdf");
 

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Diffractive transverse cross section;x;cross section (mb)");

  ofstream T_output_file("data/diff_T_sigma_x.txt");
  T_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting T integration" << endl;
  
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, T_sigma_values[i], T_sigma_errors[i]);
      T_threads[i] = thread(integrate_for_T_sigma, par);
      //this_thread::sleep_for(30s);
    }

    for (int j=0; j<x_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, T_x_values, T_sigma_values, T_x_errors, T_sigma_errors);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    T_graphs->Add(subgraph);

    for (int i=0; i<x_steps; i++) {
      ostringstream x;
      x << T_x_values[i];
      ostringstream sigma;
      sigma << T_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << T_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + x.str() + ";" + sigma.str() + ";" + sigma_err.str();
      T_output_file << line << endl;
    }
  }
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  T_sigma_canvas->Print("figures/diff_T_sigma_x_distribution.pdf");

  return 0;
  */
}
