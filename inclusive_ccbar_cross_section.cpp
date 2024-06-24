#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

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
  return normalization*L_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

double get_x(int step, int steps, double start, double stop) {
  
}

int main() {

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const int Q2_values[] = {1, 2, 3, 4, 5};

  const int x_steps = 100;
  const double x_start = 0.1;
  const double x_stop = 1e-5;
  const double x_step = log(x_stop - x_start)/x_steps;

  const bool print = false;

  const int dim = 3;
  double res, err;

  double xl[3] = {-1*integration_radius, -1*integration_radius, 0};
  double xu[3] = {integration_radius, integration_radius, 1};

  struct parameters params = {1, 1};

  if (print) {
    cout << "normalization: " << normalization << endl;
    cout << "extreme L_integrand: " << L_integrand(-100, -100, 0.5, params.Q2, params.x) << endl;
    cout << "extreme T_integrand: " << T_integrand(-100, -100, 0.5, params.Q2, params.x) << endl;
  }

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};
  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  TGraph* L_graphs[size(Q2_values)];
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];

    double L_x_values[x_steps], L_sigma_values[x_steps];
    for (int i=0; i<x_steps; i++) {
      params.x = x_start + exp(i*x_step);
      L_x_values[i] = params.x;

      cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << res << endl;

      gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);

      gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);

      if (print) {
        cout << "L warmup" << endl;
        cout << "res: " << res << endl;
        cout << "err: " << err << endl;
        cout << "chisq: " << gsl_monte_vegas_chisq(L_s) << endl;
        cout << endl;
      }

      for (int i=0; i<integration_iterations; i++) {
        gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);

        if (print) {
          cout << "res: " << res << endl;
          cout << "err: " << err << endl;
          cout << "chisq: " << gsl_monte_vegas_chisq(L_s) << endl;
          cout << endl;
        }
      }
      L_sigma_values[i] = res;

      gsl_monte_vegas_free(L_s);
    }

    L_graphs[j] = new TGraph(x_steps, L_x_values, L_sigma_values);
    L_graphs[j]->SetLineColor(j+1);
    L_graphs[j]->SetMarkerColor(j+1);
    L_graphs[j]->SetTitle("Longitudinal cross section;x ;cross section (mb)");
  }

  TCanvas* L_sigma_canvas = new TCanvas("L_sigma_canvas", "", 1000, 600);

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    if (j==0) {
      L_graphs[size(Q2_values)-1-j]->Draw("AC*");
    } else {
      L_graphs[size(Q2_values)-1-j]->Draw("C*");
    }
  }

  TLegend L_legend(.7,.6,.9,.9,"");
  L_legend.SetFillColor(0);
  L_legend.SetTextSize(0.05);

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    TString label = "Q^{2}=" + to_string(Q2_values[j]);
    L_legend.AddEntry(L_graphs[j], label);
  }
  L_legend.DrawClone("Same");

  gPad->SetLogx();

  L_sigma_canvas->Print("figures/L_sigma_x_distribution.pdf");


  TGraph* T_graphs[size(Q2_values)];

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    params.Q2 = Q2_values[j];

    double T_x_values[x_steps], T_sigma_values[x_steps];
    for (int i=0; i<x_steps; i++) {
      params.x = x_start + exp(i*x_step);
      T_x_values[i] = params.x;

      cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << res << endl;

      gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);

      gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);

      if (print) {
        cout << "T warmup" << endl;
        cout << "res: " << res << endl;
        cout << "err: " << err << endl;
        cout << "chisq: " << gsl_monte_vegas_chisq(T_s) << endl;
        cout << endl;
      }

      for (int i=0; i<integration_iterations; i++) {
        gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);

        if (print) {
          cout << "res: " << res << endl;
          cout << "err: " << err << endl;
          cout << "chisq: " << gsl_monte_vegas_chisq(T_s) << endl;
          cout << endl;
        }
      }
      T_sigma_values[i] = res;

      gsl_monte_vegas_free(T_s);
    }

    T_graphs[j] = new TGraph(x_steps, T_x_values, T_sigma_values);
    T_graphs[j]->SetLineColor(j+1);
    T_graphs[j]->SetMarkerColor(j+1);
    T_graphs[j]->SetTitle("Transverse cross section;x ;cross section (mb)");
  }

  TCanvas* T_sigma_canvas = new TCanvas("T_sigma_canvas", "", 1000, 600);

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    if (j==0) {
      T_graphs[j]->Draw("AC*");
    } else {
      T_graphs[j]->Draw("C*");
    }
  }

  TLegend T_legend(.7,.6,.9,.9,"");
  T_legend.SetFillColor(0);
  T_legend.SetTextSize(0.05);
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    TString label = "Q^{2}=" + to_string(Q2_values[j]);
    T_legend.AddEntry(T_graphs[j], label);
  }
  T_legend.DrawClone("Same");

  gPad->SetLogx();

  T_sigma_canvas->Print("figures/T_sigma_x_distribution.pdf");

  gsl_rng_free(rng);

  return 0;
}
