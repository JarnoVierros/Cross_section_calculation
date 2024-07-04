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
using namespace std;


const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double sigma_0 = 29.9416; //mb
const double B_D = sigma_0/(4*M_PI); //mb
const double Q_0 = 1; //GeV
const double x_0 = 0.000041;
const double lambda_star = 0.288;

const double normalization = 2*N_c*alpha_em*B_D*B_D/(2*M_PI*2*M_PI)*e_f*e_f;

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
  return r*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_J0(r/2*sqrt(gsl_pow_2(x0-x1)+gsl_pow_2(y0-y1)))*dipole_amplitude(r, x);
}

struct core_parameters {double x0; double y0; double x1; double y1; double z; double Q2; double x;};

double L_core_g(double r, void * params) {
  struct core_parameters *par = (struct core_parameters *)params;
  return L_core(r, par->x0, par->y0, par->x1, par->y1, par->z, par->Q2, par->x);
}

double L_integrand(double x0, double y0, double x1, double y1, double z, double Q2, double x, int seed) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;
  struct core_parameters params = {x0, y0, x1, y1, z, Q2, x};

  gsl_function F;
  F.function = &L_core_g;
  F.params = &params;

  int status = gsl_integration_qags (&F, 0, 10, 0, 0.000001, 1000, w, &result, &error);
  if (status != 0) {
    if (status == 18) {
    } else {
      cout<<status<<endl;throw (status);
    }
  }
  /*
  if (status != 0) {cout << "status: " << status << endl;
  printf ("result          = % .18f\n", result);
  printf ("estimated error = % .18f\n", error);
  printf ("intervals       = %zu\n", w->size);
  }
  */
  gsl_integration_workspace_free (w);

  //cout << x0 << "," << y0 << "," << x1 << "," << y1 << "," << z << " " << "res: " << core_res << ", err: " << core_err << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  double core_value = result;
  return 4*z*(1-z)*epsilon2(z, Q2)*exp(-B_D*(gsl_pow_2(x0+x1)+gsl_pow_2(y0+y1)))*gsl_pow_2(core_value);
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], k[4], par->Q2, par->x, 0);
}

struct thread_par_struct
{
  double Q2;
  double x;
  double &output;
  double &L_sigma_error;
  thread_par_struct(double a1, double a2, double &a3, double &a4) : Q2(a1), x(a2), output(a3), L_sigma_error(a4) {}
};

void integrate_for_sigma(thread_par_struct par) {

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 1000000;
  const int integration_iterations = 1;

  const int dim = 5;
  double res, err;

  double xl[5] = {-1*integration_radius, -1*integration_radius, -1*integration_radius, -1*integration_radius, 0};
  double xu[5] = {integration_radius, integration_radius, integration_radius, integration_radius, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  double &output = par.output;
  double &L_sigma_error = par.L_sigma_error;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
  if (status != 0) {cout<<status<<endl;throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
    if (status != 0) {cout<<status<<endl;throw (status);}
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x=" << params.x << endl;
  }
  output = res;
  L_sigma_error = err;
  cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << res << ", err: " << err << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
  if (status != 0) {cout<<status<<endl;throw (status);}
}

int main() {

  gsl_set_error_handler_off();

  if (false) {
    for (int i=0; i<1; i++) {
      L_integrand(0.1, 0.2, 0.3, 0.4, 0.001, 1, 0.001, i+1);
    }
    return 0;
  }

  const int Q2_values[] = {1};

  const int x_steps = 8;
  const double x_start = 1e-5;
  const double x_stop = 0.1;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Longitudinal cross section;x;cross section (GeV^(-2))");
  cout << "Start integration" << endl;
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double L_x_values[x_steps], L_sigma_values[x_steps], L_x_errors[x_steps], L_sigma_errors[x_steps];
    thread threads[x_steps];
    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      L_x_values[i] = x;
      L_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, L_sigma_values[i], L_sigma_errors[i]);
      threads[i] = thread(integrate_for_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      threads[j].join();
    }
    TGraphErrors* subgraph = new TGraphErrors(x_steps, L_x_values, L_sigma_values, L_x_errors, L_sigma_errors);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);
  }

  TCanvas* L_sigma_canvas = new TCanvas("diff_L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  gPad->SetLogy();

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
      if (status != 0) {cout<<status<<endl;throw (status);}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
        if (status != 0) {cout<<status<<endl;throw (status);}
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
  
  return 0;
}
