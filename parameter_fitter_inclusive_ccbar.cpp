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
#include "TMinuit.h"

#include <string>
#include <iostream>
#include <fstream>
using namespace std;


const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double Q_0 = 1; //GeV

const double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

static double sigma_0 = 29.12; //mb
static double x_0 = 0.000041;
static double lambda_star = 0.288;

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

struct parameters {double Q2; double x; double y;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

static vector<double> Q2_values;
static vector<double> x_values;
static vector<double> y_values;
static vector<double> measured_sigma_values;
void data_fit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  double chisq = 0;
  double delta;

  sigma_0 = par[0];
  x_0 = par[1];
  lambda_star = par[2];

  const int dim = 3;
  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_iterations = 1;
  const int integration_calls = 100000;
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
    double sigma_r = F_2 - y_values[j]*y_values[j]/(1+gsl_pow_2(1-y_values[j]))*F_L;
  
    delta = (sigma_r - measured_sigma_values[j])/1; //put error in denominator
    chisq += delta*delta;
  }
  f = chisq;
  gsl_rng_free(rng);
}

int main() {

  const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  ifstream data_file("HERA_data.dat");

  string line;
  while(getline (data_file, line)) {
    int i = 0;
    string value = "";
    while(line[i] != ' ') {
      value += line[i];
    }
    Q2_values.push_back(stod(value));

    value = "";
    while(line[i] != ' ') {
      value += line[i];
    }
    x_values.push_back(stod(value));

    value = "";
    while(line[i] != ' ') {
      value += line[i];
    }
    y_values.push_back(stod(value));

    value = "";
    while(line[i] != ' ') {
      value += line[i];
    }
    measured_sigma_values.push_back(stod(value));
    value = "";
  }

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  Double_t val1,err1,val2,err2,val3,err3;
  Double_t arglist[10];
  Int_t ierflg = 0;
  static Double_t vstart[7];
  static Double_t step[7];

  TMinuit* gMinuit = new TMinuit(3);
  gMinuit->SetFCN(data_fit);

  arglist[0] = 1;

  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  vstart[0] = sigma_0;
  vstart[1] = x_0;
  vstart[2] = lambda_star;
  step[0] = 1;
  step[1] = 0.000001;
  step[2] = 0.01;
  gMinuit->mnparm(0, "a0", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
  gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);

  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  gMinuit->GetParameter(0, val1, err1);
  gMinuit->GetParameter(1, val2, err2);
  gMinuit->GetParameter(2, val3, err3);

  cout << "val1: " << val1 << ", err1: " << err1 << endl;
  cout << "val2: " << val1 << ", err2: " << err1 << endl;
  cout << "val3: " << val1 << ", err3: " << err1 << endl;
  
  return 0;
}
