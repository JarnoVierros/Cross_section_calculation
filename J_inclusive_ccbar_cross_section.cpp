#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
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

int main() {

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const int Q2_values[] = {1, 3, 5, 8, 10};

  const int x_steps = 30;
  const double x_start = 1e-5;
  const double x_stop = 0.01;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  string filename = "data/dipole_amplitude_with_IP_dependence.csv";
  load_dipole_amplitudes(table, filename);

  const int dim = 4;
  double res, err;

  double xl[4] = {0, 0, 0, 0};
  double xu[4] = {34.64, 17.32, M_PI, 1};

  struct parameters params = {1, 1};

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};
  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Longitudinal cross section;x;cross section (mb)");

  ofstream L_output_file("data/J_inclusive_L_sigma_x.txt");
  L_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];
    double L_x_values[x_steps], L_sigma_values[x_steps], L_sigma_errors[x_steps];
    for (int i=0; i<x_steps; i++) {

      params.x = pow(10, log10(x_start) + i*x_step);
      L_x_values[i] = params.x;

      gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }
      L_sigma_values[i] = res; // no unit change to mb needed(?)
      L_sigma_errors[i] = err;

      cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << res << endl;

      ostringstream x;
      x << L_x_values[i];
      ostringstream sigma;
      sigma << L_sigma_values[i];
      ostringstream sigma_error;
      sigma_error << L_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + x.str() + ";" + sigma.str() + ";" + sigma_error.str();
      L_output_file << line << endl;

      gsl_monte_vegas_free(L_s);
      if (status != 0) {throw "gsl error";}
    }
    TGraph* subgraph = new TGraph(x_steps, L_x_values, L_sigma_values);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);
  }
  L_output_file.close();

  TCanvas* L_sigma_canvas = new TCanvas("L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  L_sigma_canvas->Print("figures/J_L_sigma_x_distribution.pdf");
  

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Transverse cross section;x;cross section (mb)");

  ofstream T_output_file("data/J_inclusive_T_sigma_x.txt");
  T_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];

    double T_x_values[x_steps], T_sigma_values[x_steps], T_sigma_errors[x_steps];
    for (int i=0; i<x_steps; i++) {
      params.x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = params.x;

      gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }

      T_sigma_values[i] = res;
      T_sigma_errors[i] = err;

      cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << res << endl;

      ostringstream x;
      x << T_x_values[i];
      ostringstream sigma;
      sigma << T_sigma_values[i];
      ostringstream sigma_error;
      sigma_error << T_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + x.str() + ";" + sigma.str() + ";" + sigma_error.str();
      T_output_file << line << endl;

      gsl_monte_vegas_free(T_s);
    }
    TGraph* subgraph = new TGraph(x_steps, T_x_values, T_sigma_values);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    T_graphs->Add(subgraph);
  }
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  T_sigma_canvas->Print("figures/J_T_sigma_x_distribution.pdf");

  gsl_rng_free(rng);
  
  return 0;
}