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
using namespace std;


const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double sigma_0 = 2.99416e+01; //mb
const double Q_0 = 1; //GeV
const double x_0 = 7.67079e-05;
const double lambda_star = 3.64361e-01;

const double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

double epsilon(double z, double Q2) {
  return sqrt(m_f*m_f + z*(1-z)*Q2);
}

double dipole_amplitude(double r, double Q2, double W) {
  return sigma_0*(1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow((Q2)/(W*W+Q2)/x_0, lambda_star/2)))));
}

double L_integrand(double r_x, double r_y, double z, double Q2, double W) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return 4*Q2*z*z*(1-z)*(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, Q2, W);
}

double T_integrand(double r_x, double r_y, double z, double Q2, double W) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return (m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + gsl_pow_2(epsilon(z, Q2))*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, Q2, W);
}

struct parameters {double Q2; double W;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], par->Q2, par->W);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], par->Q2, par->W);
}

int main() {

  gsl_set_error_handler_off();

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const double Q2_values[] = {1e2, 1e3, 1e4, 1e5, 1e6, 1e7};

  const int W_steps = 100;
  const double W_start = 1;
  const double W_stop = 1e5;
  const double W_step = 1.0/(W_steps-1)*log10(W_stop/W_start);

  const int dim = 3;
  double res, err;

  double xl[3] = {-1*integration_radius, -1*integration_radius, 0};
  double xu[3] = {integration_radius, integration_radius, 1};

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
  L_graphs->SetTitle("Longitudinal cross section;W (GeV);cross section (mb)");
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];
    double L_W_values[W_steps], L_sigma_values[W_steps];
    for (int i=0; i<W_steps; i++) {

      params.W = pow(10, log10(W_start) + i*W_step);
      L_W_values[i] = params.W;

      cout << "L, Q²=" << params.Q2 << ", W=" << params.W << ", res: " << res << endl;

      gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }
      L_sigma_values[i] = res;

      gsl_monte_vegas_free(L_s);
      if (status != 0) {throw "gsl error";}
    }
    TGraph* subgraph = new TGraph(W_steps, L_W_values, L_sigma_values);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);
  }

  TCanvas* L_sigma_canvas = new TCanvas("L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.15, 0.55, 0.3, 0.9);

  L_sigma_canvas->Print("figures/L_sigma_W_distribution.pdf");
  

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Transverse cross section;W (GeV);cross section (mb)");

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];

    double T_W_values[W_steps], T_sigma_values[W_steps];
    for (int i=0; i<W_steps; i++) {
      params.W = pow(10, log10(W_start) + i*W_step);
      T_W_values[i] = params.W;

      cout << "T, Q²=" << params.Q2 << ", W=" << params.W << ", res: " << res << endl;

      gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }

      T_sigma_values[i] = res;

      gsl_monte_vegas_free(T_s);
    }
    TGraph* subgraph = new TGraph(W_steps, T_W_values, T_sigma_values);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    T_graphs->Add(subgraph);
  }

  TCanvas* T_sigma_canvas = new TCanvas("T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.15, 0.55, 0.3, 0.9);

  T_sigma_canvas->Print("figures/T_sigma_W_distribution.pdf");

  gsl_rng_free(rng);
  
  return 0;
}
