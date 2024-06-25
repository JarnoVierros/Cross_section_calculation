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

struct parameters {double Q2; double x; double eta;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

int main() {

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  double Q2_values[] = {1, 2, 3};
  double x_values[] = {1, 2, 3};
  double eta_values[] = {1, 2, 3};

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

  for (int j=0; j<size(Q2_values); j++) {

    params.Q2 = Q2_values[j];
    params.x = x_values[j];

    gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);

    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
    if (status != 0) {throw "gsl error";}

    for (int i=0; i<integration_iterations; i++) {
      status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
      if (status != 0) {throw "gsl error";}
    }
    
    double sigma_L = res;

    gsl_monte_vegas_free(L_s);

    gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);

    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
    if (status != 0) {throw "gsl error";}

    for (int i=0; i<integration_iterations; i++) {
      status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
      if (status != 0) {throw "gsl error";}
    }

    double sigma_T = res;

    gsl_monte_vegas_free(T_s);

    double F_L = Q2_values[j]/(4*M_PI*M_PI*alpha_em)*sigma_L;
    double F_T = Q2_values[j]/(4*M_PI*M_PI*alpha_em)*sigma_T;
    double F_2 = F_L + F_T;
    double sigma_r = F_2 - eta_values[j]*eta_values[j]/(1+gsl_pow_2(1-eta_values[j]))*F_L;

    cout << "sigma_r: " << sigma_r << endl;

  }

  gsl_rng_free(rng);
  
  return 0;
}
