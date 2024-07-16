#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <thread>
using namespace std;

#include "direct_dipole_amp_reader.h"

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double normalization = 8/M_PI*alpha_em*N_c*e_f*e_f;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double x) {
  return get_dipole_amplitude(table, r, b_min, phi, x);
}

double L_integrand(double r, double b_min, double phi, double z, double Q2, double x) {
  return r*b_min*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, b_min, phi, x);
}

double T_integrand(double r, double b_min, double phi, double z, double Q2, double x) {
  return r*b_min*(m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, b_min, phi, x);
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], par->Q2, par->x);
}

struct thread_par_struct
{
  double Q2;
  double x;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double a2, double &a3, double &a4) : Q2(a1), x(a2), sigma(a3), sigma_error(a4) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  //const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const int dim = 4;
  double res, err;

  double xl[4] = {0, 0, 0, 0};
  double xu[4] = {34.64, 17.32, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
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
    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
    if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {

  //const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const int dim = 4;
  double res, err;

  double xl[4] = {0, 0, 0, 0};
  double xu[4] = {34.64, 17.32, M_PI, 1};

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

int main() {

  gsl_set_error_handler_off();

  const int Q2_values[] = {1, 3, 5, 8, 10};

  const int x_steps = 30;
  const double x_start = 1e-5;
  const double x_stop = 0.01;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  string filename = "data/dipole_amplitude_with_IP_dependence.csv";
  load_dipole_amplitudes(table, filename);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Longitudinal cross section;x;cross section (mb)");

  ofstream L_output_file("data/J_L_inclusive_sigma_x.txt");
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

  TCanvas* L_sigma_canvas = new TCanvas("J_L_inclusive_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  L_sigma_canvas->Print("figures/J_L_inclusive_sigma_x_distribution.pdf");
 

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Transverse cross section;x;cross section (mb)");

  ofstream T_output_file("data/J_T_inclusive_sigma_x.txt");
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

  TCanvas* T_sigma_canvas = new TCanvas("J_T_inclusive_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  T_sigma_canvas->Print("figures/J_T_inclusive_sigma_x_distribution.pdf");

  return 0;
}
