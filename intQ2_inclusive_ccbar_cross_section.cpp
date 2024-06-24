#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

#include <string>
#include <iostream>
using namespace std;


const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double sigma_0 = 29.12; //mb
const double Q_0 = 1; //GeV
const double x_0 = 0.000041;
const double lambda_star = 0.288;

const double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

double epsilon(double z, double Q2) {
  return sqrt(m_f*m_f + z*(1-z)*Q2);
}

double dipole_amplitude(double r, double x) {
  return sigma_0*(1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow(x/x_0, lambda_star/2)))));
}

double L_integrand(double r_x, double r_y, double z, double Q2, double x) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return 4*Q2*z*z*(1-z)*(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, x);
}

double T_integrand(double r_x, double r_y, double z, double Q2, double x) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return (m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + gsl_pow_2(epsilon(z, Q2))*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, x);
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], par->x);
}

int main() {

  gsl_set_error_handler_off();
  int status = 0;

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const int x_steps = 10;
  const double x_start = 0.001;
  const double x_stop = 0.1;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  const int dim = 4;
  double res, err;

  double xl[4] = {-1*integration_radius, -1*integration_radius, 0, 0};
  double xu[4] = {integration_radius, integration_radius, 1, 100};

  struct parameters params = {1};

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};
  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  double L_x_values[x_steps], L_sigma_values[x_steps];
  for (int i=0; i<x_steps; i++) {

    params.x = pow(10, log10(x_start) + i*x_step);
    L_x_values[i] = params.x;

    cout << "x=" << params.x << ", res: " << res << ", integrand at max Q2: " << L_integrand(0.1, 0.1, 0.5, xu[3], params.x) << endl;

    gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);

    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
    if (status != 0) {throw "gsl error";}

    for (int i=0; i<integration_iterations; i++) {
      status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
      if (status != 0) {throw "gsl error";}
    }
    L_sigma_values[i] = res;

    gsl_monte_vegas_free(L_s);
  }
  TGraph* L_sigma_graph = new TGraph(x_steps, L_x_values, L_sigma_values);
  L_sigma_graph->SetTitle("Longitudinal cross section integrated over Q^{2}");
  L_sigma_graph->GetXaxis()->SetTitle("x");
  L_sigma_graph->GetYaxis()->SetTitle("#sigma (mb)");

  TCanvas* L_sigma_canvas = new TCanvas("L_sigma_canvas", "", 1000, 600);
  L_sigma_graph->Draw("AC*");

  gPad->SetLogx();

  L_sigma_canvas->Print("figures/L_sigma_x_distribution_integrated_Q2.pdf");
  

  double T_x_values[x_steps], T_sigma_values[x_steps];
  for (int i=0; i<x_steps; i++) {
    params.x = pow(10, log10(x_start) + i*x_step);
    T_x_values[i] = params.x;

    cout << "x=" << params.x << ", res: " << res << endl;

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
  TGraph* T_sigma_graph = new TGraph(x_steps, T_x_values, T_sigma_values);
  T_sigma_graph->SetTitle("Transverse cross section integrated over Q^{2}");

  TCanvas* T_sigma_canvas = new TCanvas("T_sigma_canvas", "", 1000, 600);
  T_sigma_graph->Draw("AC*");

  gPad->SetLogx();

  T_sigma_canvas->Print("figures/T_sigma_x_distribution_integrated_Q2.pdf");

  gsl_rng_free(rng);
  
  return 0;
}
