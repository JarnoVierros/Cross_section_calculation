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

const double normalization = 8/(2*M_PI)*alpha_em*N_c*e_f*e_f;

const double r_limit; // 34.64
const double b_min_limit; // 17.32

const int warmup_calls = 10000;
const int integration_calls = 100000;
const int integration_iterations = 1;

const string dipole_amp_type = "bfkl";
const string nucleus_type = "p";
const string filename_end = "";

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double W, double Q2) {
  return get_dipole_amplitude(table, r, b_min, phi, Q2/(W*W+Q2));
}

double L_integrand(double r, double b_min, double phi, double z, double Q2, double W) {
  return r*b_min*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*gsl_pow_2(dipole_amplitude(r, b_min, phi, W, Q2));;
}

double T_integrand(double r, double b_min, double phi, double z, double Q2, double W) {
  return r*b_min*(m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*gsl_pow_2(dipole_amplitude(r, b_min, phi, W, Q2));
}

struct parameters {double Q2; double W;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], par->Q2, par->W);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], par->Q2, par->W);
}

struct thread_par_struct
{
  double Q2;
  double W;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double a2, double &a3, double &a4) : Q2(a1), W(a2), sigma(a3), sigma_error(a4) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  const int dim = 4;
  double res, err;

  double xl[4] = {0, 0, 0, 0};
  double xu[4] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.W = par.W;
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
    cout << "nan found at x=" << params.W << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "L, Q²=" << params.Q2 << ", x=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 4;
  double res, err;

  double xl[4] = {0, 0, 0, 0};
  double xu[4] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.W = par.W;
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
    cout << "nan found at x=" << params.W << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "T, Q²=" << params.Q2 << ", x=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  if (nucleus_type == "Pb") {
    const double r_limit = 657; // 34.64
    const double b_min_limit = 328; // 17.32
  } else if (nucleus_type == "p") {
    const double r_limit = 34.64;
    const double b_min_limit = 17.32;
  } else {
    cout << "invalid nucleus type" << endl;
    throw 1;
  }

  const int Q2_values[] = {1, 3, 5, 8, 10};

  const int W_steps = 30;
  const double W_start = 1e-5;
  const double W_stop = 0.01;
  const double W_step = 1.0/(W_steps-1)*log10(W_stop/W_start);

  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
  load_dipole_amplitudes(table, filename);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Diffractive longitudinal cross section;x;cross section (mb)");

  ofstream L_output_file("data/diff_L_sigma_W_"+dipole_amp_type+"_"+nucleus_type+".txt");
  L_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting L integration" << endl;
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double L_W_values[W_steps], L_sigma_values[W_steps], L_W_errors[W_steps], L_sigma_errors[W_steps];
    thread L_threads[W_steps];

    for (int i=0; i<W_steps; i++) {
      double W = pow(10, log10(W_start) + i*W_step);
      L_W_values[i] = W;
      L_W_errors[i] = 0;
      thread_par_struct par(Q2_values[j], W, L_sigma_values[i], L_sigma_errors[i]);
      L_threads[i] = thread(integrate_for_L_sigma, par);
    }

    for (int j=0; j<W_steps; j++) {
      L_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(W_steps, L_W_values, L_sigma_values, L_W_errors, L_sigma_errors);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);

    for (int i=0; i<W_steps; i++) {
      ostringstream W;
      W << L_W_values[i];
      ostringstream sigma;
      sigma << L_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << L_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + W.str() + ";" + sigma.str() + ";" + sigma_err.str();
      L_output_file << line << endl;
    }
  }
  L_output_file.close();

  TCanvas* L_sigma_canvas = new TCanvas("diff_L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  TString fig_filename = "figures/diff_L_sigma_W_"+dipole_amp_type+"_"+nucleus_type+".pdf";
  L_sigma_canvas->Print(fig_filename);
 

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Diffractive transverse cross section;x;cross section (mb)");

  ofstream T_output_file("data/diff_T_sigma_W_"+dipole_amp_type+"_"+nucleus_type+".txt");
  T_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting T integration" << endl;
  
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_W_values[W_steps], T_sigma_values[W_steps], T_W_errors[W_steps], T_sigma_errors[W_steps];
    thread T_threads[W_steps];

    for (int i=0; i<W_steps; i++) {
      double W = pow(10, log10(W_start) + i*W_step);
      T_W_values[i] = W;
      T_W_errors[i] = 0;
      thread_par_struct par(Q2_values[j], W, T_sigma_values[i], T_sigma_errors[i]);
      T_threads[i] = thread(integrate_for_T_sigma, par);
      //this_thread::sleep_for(30s);
    }

    for (int j=0; j<W_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(W_steps, T_W_values, T_sigma_values, T_W_errors, T_sigma_errors);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    T_graphs->Add(subgraph);

    for (int i=0; i<W_steps; i++) {
      ostringstream W;
      W << T_W_values[i];
      ostringstream sigma;
      sigma << T_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << T_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + W.str() + ";" + sigma.str() + ";" + sigma_err.str();
      T_output_file << line << endl;
    }
  }
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  fig_filename = "figures/diff_T_sigma_W_"+dipole_amp_type+"_"+nucleus_type+".pdf";
  T_sigma_canvas->Print(fig_filename);

  return 0;
}
