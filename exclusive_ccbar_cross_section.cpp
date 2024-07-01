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

const double B_D = 6; //GeV^(-2)
const double Q_0 = 1; //GeV
const double x_0 = 0.000041;
const double lambda_star = 0.288;

const double normalization = 8*N_c*alpha_em/(2*M_PI*2*M_PI)*e_f*e_f*B_D;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double x) {
  return (1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow(x/x_0, lambda_star/2)))));
}

double L_core(double r, double x0, double y0, double x1, double y1, double z, double Q2, double x) {
  return r*gsl_sf_bessel_J0(r/2*sqrt(gsl_pow_2(x0-x1)+gsl_pow_2(y0-y1)))*gsl_sf_bessel_K0(epsilon(x, Q2)*r)*dipole_amplitude(r, x);
}

struct core_parameters {double x0; double y0; double x1; double y1; double z; double Q2; double x;};

double L_core_g(double *k, size_t dim, void * params) {
  struct core_parameters *par = (struct core_parameters *)params;
  return L_core(k[0], par->x0, par->y0, par->x1, par->y1, par->z, par->Q2, par->x);
}

double L_integrand(double x0, double y0, double x1, double y1, double z, double Q2, double x, int seed) {

  const double core_integration_radius = 10;
  const int core_warmup_calls = 100;
  const int core_integration_calls = 3000;
  const int core_integration_iterations = 1;
  const int core_dim = 1;

  double core_res, core_err;

  double xl[1] = {0};
  double xu[1] = {core_integration_radius};

  struct core_parameters params = {x0, y0, x1, y1, z, Q2, x};

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_core_G = {&L_core_g, core_dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  //gsl_rng_set(rng, seed);

  gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(core_dim);

  status = gsl_monte_vegas_integrate(&L_core_G, xl, xu, core_dim, core_warmup_calls, rng, L_s, &core_res, &core_err);
  if (status != 0) {throw "gsl error";}

  for (int i=0; i<core_integration_iterations; i++) {
    status = gsl_monte_vegas_integrate(&L_core_G, xl, xu, core_dim, core_integration_calls, rng, L_s, &core_res, &core_err);
    if (status != 0) {throw "gsl error";}
  }

  //cout << x0 << "," << y0 << "," << x1 << "," << y1 << "," << z << " " << "res: " << core_res << ", err: " << core_err << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  double core_value = core_res;

  gsl_monte_vegas_free(L_s);
  gsl_rng_free(rng);
  
  return z*(1-z)*epsilon2(z, Q2)*exp(-(gsl_pow_2(x0+x1)+gsl_pow_2(y0+y1))*B_D)*gsl_pow_2(core_value);
}

struct parameters {double Q2; double x;};

int integration_index = 0;
double L_g(double *k, size_t dim, void * params) {
  if (integration_index%1000==0){
    cout << integration_index << endl;
  }
  integration_index++;
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], k[4], par->Q2, par->x, 0);
}

int main() {

  gsl_set_error_handler_off();

  if (false) {
    for (int i=0; i<10; i++) {
      L_integrand(90, 91, 92, 93, 0.001, 1, 0.001, i+1);
    }
    return 0;
  }
  
  
  const double integration_radius = 100;
  const int warmup_calls = 100;
  const int integration_calls = 10000;
  const int integration_iterations = 1;

  const int Q2_values[] = {1};

  const int x_steps = 10;
  const double x_start = 1e-5;
  const double x_stop = 0.1;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  const int dim = 5;
  double res, err;

  double xl[5] = {-1*integration_radius, -1*integration_radius, -1*integration_radius, -1*integration_radius, 0};
  double xu[5] = {integration_radius, integration_radius, integration_radius, integration_radius, 1};

  struct parameters params = {1, 1};

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Longitudinal cross section;x;cross section (GeV^(-2))");
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];
    double L_x_values[x_steps], L_sigma_values[x_steps];
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
      if (gsl_isnan(res)) {
        res = 0;
        cout << "nan found at x=" << params.x << endl;
      }
      L_sigma_values[i] = res;

      cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << res << ", err: " << err << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;
      integration_index = 0;

      gsl_monte_vegas_free(L_s);
      if (status != 0) {throw "gsl error";}
    }
    TGraph* subgraph = new TGraph(x_steps, L_x_values, L_sigma_values);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);
  }

  TCanvas* L_sigma_canvas = new TCanvas("diff_L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  L_sigma_canvas->Print("figures/diff_L_sigma_x_distribution.pdf");
  
  /*
  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Transverse cross section;x;cross section (GeV^(-2))");

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];

    double T_x_values[x_steps], T_sigma_values[x_steps];
    for (int i=0; i<x_steps; i++) {
      params.x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = params.x;

      cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << res << endl;

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
    TGraph* subgraph = new TGraph(x_steps, T_x_values, T_sigma_values);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    T_graphs->Add(subgraph);
  }

  TCanvas* T_sigma_canvas = new TCanvas("T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  T_sigma_canvas->Print("figures/T_sigma_x_distribution.pdf");

  */

  gsl_rng_free(rng);
  
  return 0;
}
