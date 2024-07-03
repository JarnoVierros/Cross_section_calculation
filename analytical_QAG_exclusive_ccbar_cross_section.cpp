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
using namespace std;


const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double sigma_0 = 29.9416; //GeV^(-2)
const double Q_0 = 1; //GeV
const double x_0 = 7.67079e-05;
const double lambda_star = 3.64361e-01;

const double normalization = N_c*alpha_em*e_f*e_f*sigma_0/(2*gsl_pow_2(2*M_PI));

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double x) {
  return (1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow(x/x_0, lambda_star/2)))));
}

double L_core(double r, double px, double py, double z, double Q2, double x) {
  return r*gsl_sf_bessel_J0(r/2*sqrt(px*px+py*py))*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*dipole_amplitude(r, x);
}

double T_core_1(double r, double px, double py, double z, double Q2, double x) {
  return r*epsilon(z, Q2)*gsl_sf_bessel_J1(r/2*sqrt(px*px+py*py))*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*dipole_amplitude(r, x);
}

double T_core_2(double r, double px, double py, double z, double Q2, double x) {
  return r*gsl_sf_bessel_J0(r/2*sqrt(px*px+py*py))*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*dipole_amplitude(r, x);
}

struct core_parameters {double px; double py; double z; double Q2; double x;};

double L_core_g(double r, void * params) {
  struct core_parameters *par = (struct core_parameters *)params;
  return L_core(r, par->px, par->py, par->z, par->Q2, par->x);
}

double T_core_1_g(double r, void * params) {
  struct core_parameters *par = (struct core_parameters *)params;
  return T_core_1(r, par->px, par->py, par->z, par->Q2, par->x);
}

double T_core_2_g(double r, void * params) {
  struct core_parameters *par = (struct core_parameters *)params;
  return T_core_2(r, par->px, par->py, par->z, par->Q2, par->x);
}

double L_integrand(double px, double py, double z, double Q2, double x, int seed) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;
  struct core_parameters params = {px, py, z, Q2, x};

  gsl_function F;
  F.function = &L_core_g;
  F.params = &params;

  int status = gsl_integration_qags (&F, 0, 10, 0, 0.001, 1000, w, &result, &error);
  if (status != 0) {
    if (status == 18) {
    } else if (status == 22) {
      cout << "warning: divergent qags integral" << endl;
    } else {
      cout << "L_integrand: " << status << endl;
      cout << px << "," << py << "," << z << "," << Q2 << "," << x << " " << "res: " << result << ", err: " << error << endl;
      throw (status);
    }
  }
  /*
  printf ("result          = % .18f\n", result);
  printf ("estimated error = % .18f\n", error);
  printf ("intervals       = %zu\n", w->size);
  */
  gsl_integration_workspace_free (w);

  //cout << px << "," << py << "," << z << "," << Q2 << "," << x << " " << "res: " << result << ", err: " << error << endl;

  double core_value = result;
  return 4*z*(1-z)*epsilon2(z, Q2)*gsl_pow_2(core_value);
}

struct core_integration_thread_struct
{
  double px;
  double py;
  double z;
  double Q2;
  double x;
  double &output;
  core_integration_thread_struct(double a1, double a2, double a3, double a4, double a5, double &a6) : px(a1), py(a2), z(a3), Q2(a4), x(a5), output(a6) {}
};

void integrate_T_core_1(core_integration_thread_struct par) {
  //this_thread::sleep_for(10s);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;
  struct core_parameters params = {par.px, par.py, par.z, par.Q2, par.x};
  double &output = par.output;

  gsl_function F;
  F.function = &T_core_1_g;
  F.params = &params;
  int status = gsl_integration_qags (&F, 0, 10, 0, 0.001, 1000, w, &result, &error);
  if (status != 0) {
    if (status == 18) {
    } else if (status == 22) {
      cout << "warning: divergent qags integral" << endl;
    } else {
      cout << "integrate_T_core_1 status is: " << status << endl;
      throw (status);
    }
  }
  /*
  printf ("result          = % .18f\n", result);
  printf ("estimated error = % .18f\n", error);
  printf ("intervals       = %zu\n", w->size);
  */
  
  gsl_integration_workspace_free (w);

  //cout << par.px << "," << par.py << "," << par.z << "," << par.Q2 << "," << par.x << " " << "res: " << result << ", err: " << error << endl;

  output = result;
}

void integrate_T_core_2(core_integration_thread_struct par) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error;
  struct core_parameters params = {par.px, par.py, par.z, par.Q2, par.x};
  double &output = par.output;

  gsl_function F;
  F.function = &T_core_2_g;
  F.params = &params;

  int status = gsl_integration_qags (&F, 0, 10, 0, 0.001, 1000, w, &result, &error);
  if (status != 0) {
    if (status == 18) {
    } else if (status == 22) {
      cout << "warning: divergent qags integral" << endl;
    } else {
      cout << "integrate_T_core_2 exception: " << status << endl;
      throw (status);
    }
  }
  /*
  cout << "integrate_T_core_2 results" << endl;
  printf ("result          = % .18f\n", result);
  printf ("estimated error = % .18f\n", error);
  printf ("intervals       = %zu\n", w->size);
  */
  gsl_integration_workspace_free (w);

  //cout << par.px << "," << par.py << "," << par.z << "," << par.Q2 << "," << par.x << " " << "res: " << result << ", err: " << error << endl;

  output = result;
}

double T_integrand(double px, double py, double z, double Q2, double x, int seed) {
  double core_1_integral;
  core_integration_thread_struct core_1_parameters(px, py, z, Q2, x, core_1_integral);
  //thread core_1_thread(integrate_T_core_1, core_1_parameters);
  integrate_T_core_1(core_1_parameters);
  double core_2_integral;
  core_integration_thread_struct core_2_parameters(px, py, z, Q2, x, core_2_integral);
  //thread core_2_thread(integrate_T_core_2, core_2_parameters);
  integrate_T_core_2(core_2_parameters);
  //core_1_thread.join();
  //core_2_thread.join();

  return (z*z + gsl_pow_2(1-z))*gsl_pow_2(core_1_integral) + m_f*m_f*gsl_pow_2(core_2_integral);
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], par->Q2, par->x, 0);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], par->Q2, par->x, 0);
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

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const int dim = 3;
  double res, err;

  double xl[3] = {0, 0, 0};
  double xu[3] = {integration_radius, integration_radius, 1};

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

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const int dim = 3;
  double res, err;

  double xl[3] = {0, 0, 0};
  double xu[3] = {integration_radius, integration_radius, 1};

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

  if (false) {
    for (int i=0; i<1; i++) {
      L_integrand(0.1, 0.2, 0.001, 1, 0.001, i+1);
    }
    return 0;
  }

  const int Q2_values[] = {1, 3, 5, 8, 10};

  const int x_steps = 10;
  const double x_start = 1e-5;
  const double x_stop = 0.1;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Diffractive longitudinal cross section;x;cross section (mb)");
  cout << "Starting L integration" << endl;
  double L_x_values[x_steps], L_sigma_values[x_steps], L_x_errors[x_steps], L_sigma_errors[x_steps];
  for (long unsigned int j=0; j<size(Q2_values); j++) {
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
  }

  TCanvas* L_sigma_canvas = new TCanvas("diff_L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  L_sigma_canvas->Print("figures/analytical_diff_L_sigma_x_distribution.pdf");


  ofstream L_output_file("diff_L_sigma_x.txt");

  L_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  for (long unsigned int i=0; i<size(Q2_values); i++) {
    for (int j=0; j<x_steps; j++) {
      string line = to_string(Q2_values[i]) + ";" + to_string(L_x_values[j]) + ";" + to_string(L_sigma_values[j]) + ";" + to_string(L_sigma_errors[j]);
      L_output_file << line << endl;
    }
  }
  L_output_file.close();

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Diffractive transverse cross section;x;cross section (mb)");
  cout << "Starting T integration" << endl;
  double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps];
  for (long unsigned int j=0; j<size(Q2_values); j++) {
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
  }

  TCanvas* T_sigma_canvas = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  T_sigma_canvas->Print("figures/analytical_diff_T_sigma_x_distribution.pdf");
  
  
  ofstream T_output_file("diff_T_sigma_x.txt");

  T_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  for (long unsigned int i=0; i<size(Q2_values); i++) {
    for (int j=0; j<x_steps; j++) {
      string line = to_string(Q2_values[i]) + ";" + to_string(T_x_values[j]) + ";" + to_string(T_sigma_values[j]) + ";" + to_string(T_sigma_errors[j]);
      L_output_file << line << endl;
    }
  }
  T_output_file.close();

  return 0;
}
