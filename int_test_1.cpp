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
const double m = 1.27; //GeV

const double sigma_0 = 29.12; //mb
const double Q_0 = 1; //GeV
const double x_0 = 0.000041;
const double lambda_star = 0.288;

double Q;
double x;

double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

double epsilon(double z) {
  return sqrt(m*m + z*(1-z)*Q*Q);
}

double dipole_amplitude(double r) {
  return sigma_0*(1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow(x/x_0, lambda_star/2)))));
}

double L_integrand(double r_x, double r_y, double z) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return 4*Q*Q*z*z*(1-z)*(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z)*r))*dipole_amplitude(r);
}

double T_integrand(double r_x, double r_y, double z) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return (m*m*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z)*r)) + gsl_pow_2(epsilon(z))*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z)*r)))*dipole_amplitude(r);
}

double L_g(double *k, size_t dim, void *params) {
  return normalization*L_integrand(k[0], k[1], k[2]);
}

double T_g(double *k, size_t dim, void *params) {
  return normalization*T_integrand(k[0], k[1], k[2]);
}

int main() {

  bool print = false;
  int iterations = 1;

  const int dim = 3;
  double res, err;

  double xl[3] = {-100, -100, 0};
  double xu[3] = {100, 100, 1};

  if (print) {
    cout << "normalization: " << normalization << endl;
    cout << "extreme L_integrand: " << L_integrand(-100, -100, 0.5) << endl;
    cout << "extreme T_integrand: " << T_integrand(-100, -100, 0.5) << endl;
  }

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, 0};
  gsl_monte_function T_G = {&T_g, dim, 0};

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  TGraph* L_grapsh[5];
  for (int j=0; j<5; j++) {
    Q = 1 + j;

    double L_x_values[100], L_sigma_values[100];
    for (int i=0; i<100; i++) {
      x = (i+1)*0.001;
      L_x_values[i] = x;

      cout << "L, Q=" << Q << ", x=" << x << ", res: " << res << endl;

      gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);

      gsl_monte_vegas_integrate(&L_G, xl, xu, dim, 10000, rng, L_s, &res, &err);

      if (print) {
        cout << "L warmup" << endl;
        cout << "res: " << res << endl;
        cout << "err: " << err << endl;
        cout << "chisq: " << gsl_monte_vegas_chisq(L_s) << endl;
        cout << endl;
      }

      for (int i=0; i<iterations; i++) {
        gsl_monte_vegas_integrate(&L_G, xl, xu, dim, 100000, rng, L_s, &res, &err);

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

    L_graphs[j] = new TGraph(100, L_x_values, L_sigma_values);
    L_graphs[j]->SetLineColor(j+1);
    L_graphs[j]->SetMarkerColor(j+1);
    L_graphs[j]->SetTitle("Longitudinal cross section;x;sigma");
  }

  TCanvas* L_sigma_canvas = new TCanvas("L_sigma_canvas", "", 1000, 600);

  for (int j=0; j<5; j++) {
    if (j==0) {
      L_graphs[j]->Draw("AC*");
    } else {
      L_graphs[j]->Draw("C*");
    }
  }

  TLegend L_legend(.7,.6,.9,.9,"");
  L_legend.SetFillColor(0);
  L_legend.SetTextSize(0.05);
  L_legend.AddEntry(L_graphs[0],"Q=1");
  L_legend.AddEntry(L_graphs[1],"Q=2");
  L_legend.AddEntry(L_graphs[2],"Q=3");
  L_legend.AddEntry(L_graphs[3],"Q=4");
  L_legend.AddEntry(L_graphs[4],"Q=5");
  leg1.DrawClone("Same");

  L_sigma_canvas->Print("L_sigma_x_distribution.pdf");


  TGraph* T_graphs[5];

  for (int j=0; j<5; j++) {
    Q = 1 + j;

    double T_x_values[100], T_sigma_values[100];
    for (int i=0; i<100; i++) {
      x = (i+1)*0.001;
      T_x_values[i] = x;

      cout << "T, Q=" << Q << ", x=" << x << ", res: " << res << endl;

      gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);

      gsl_monte_vegas_integrate(&T_G, xl, xu, dim, 10000, rng, T_s, &res, &err);

      if (print) {
        cout << "T warmup" << endl;
        cout << "res: " << res << endl;
        cout << "err: " << err << endl;
        cout << "chisq: " << gsl_monte_vegas_chisq(T_s) << endl;
        cout << endl;
      }

      for (int i=0; i<iterations; i++) {
        gsl_monte_vegas_integrate(&T_G, xl, xu, dim, 100000, rng, T_s, &res, &err);

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

    T_graphs[j] = new TGraph(100, T_x_values, T_sigma_values);
    T_graphs[j]->SetLineColor(j+1);
    T_graphs[j]->SetMarkerColor(j+1);
    T_graphs[j]->SetTitle("Transverse cross section;x;sigma");
  }

  T_legend.SetFillColor(0);
  T_legend.SetTextSize(0.05);
  T_legend.AddEntry(T_graphs[0],"Q=1");
  T_legend.AddEntry(T_graphs[1],"Q=2");
  T_legend.AddEntry(T_graphs[2],"Q=3");
  T_legend.AddEntry(T_graphs[3],"Q=4");
  T_legend.AddEntry(T_graphs[4],"Q=5");
  leg1.DrawClone("Same");

  TCanvas* T_sigma_canvas = new TCanvas("T_sigma_canvas", "", 1000, 600);

  for (int j=0; j<5; j++) {
    if (j==0) {
      T_graphs[5-j]->Draw("AC*");
    } else {
      T_graphs[5-j]->Draw("C*");
    }
  }

  T_sigma_canvas->Print("T_sigma_x_distribution.pdf");

  gsl_rng_free(rng);

  return 0;
}
