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
#include <thread>
using namespace std;


const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double Q_0 = 1; //GeV

const double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

static double sigma_0 = 2.99416e+01; //mb
static double x_0 = 7.67079e-05;
static double lambda_star = 3.64361e-01;

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

struct par_struct
{
  double Q2;
  double x;
  double y;
  double measured_sigma;
  double relative_measurement_error;
  double &output;
  par_struct(double a1, double a2, double a3, double a4, double a5, double &a6) : Q2(a1), x(a2), y(a3), measured_sigma(a4), relative_measurement_error(a5), output(a6) {}
};

void integrate_for_delta(par_struct par) {
  double Q2 = par.Q2;
  double x = par.x;
  double y = par.y;
  double measured_sigma = par.measured_sigma;
  double relative_measurement_error = par.relative_measurement_error;
  double &output = par.output;

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

  params.Q2 = Q2;
  params.x = x;

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
  gsl_rng_free(rng);

  double F_L = Q2/(4*M_PI*M_PI*alpha_em)*sigma_L;
  double F_T = Q2/(4*M_PI*M_PI*alpha_em)*sigma_T;
  double F_2 = F_L + F_T;
  double sigma_r = F_2 - y*y/(1+gsl_pow_2(1-y))*F_L;

  double delta = (sigma_r - measured_sigma)/(relative_measurement_error/100*measured_sigma);

  output = delta;
}

static vector<double> Q2_values;
static vector<double> x_values;
static vector<double> y_values;
static vector<double> measured_sigma_values;
static vector<double> relative_measurement_errors;
void data_fit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  double chisq = 0;

  sigma_0 = par[0];
  x_0 = par[1];
  lambda_star = par[2];

  double deltas[size(Q2_values)];
  thread threads[size(Q2_values)];
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    par_struct par(Q2_values[j], x_values[j], y_values[j], measured_sigma_values[j], relative_measurement_errors[j], deltas[j]);
    threads[j] = thread(integrate_for_delta, par);
  }

  for (long unsigned int j=0; j<size(Q2_values); j++) {
    threads[j].join();
    double delta = deltas[j];
    chisq += delta*delta;
  }

  f = chisq;
  
  cout << "sigma_0: " << sigma_0 << ", x_0: " << x_0 << ", lambda_star: " << lambda_star << ", chisq: " << chisq << endl;
}

int main() {

  gsl_set_error_handler_off();

  const string filename = "data/HERA_data.dat";

  ifstream data_file(filename);

  cout << "Reading: " << filename << endl;
  string line;
  while(getline (data_file, line)) {
    long unsigned int i = 0;
    string value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    Q2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    x_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    y_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    measured_sigma_values.push_back(stod(value));
    i++;

    for (int j=0; j<3; j++) {
      while(line[i] != ' ') {
        i++;
      }
      i++;
    }

    value = "";
    while(i<line.length()) {
      value += line[i];
      i++;
    }
    relative_measurement_errors.push_back(stod(value));
  }
  cout << "Finished reading file" << endl;

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

  vstart[0] = 2.99415e+01;
  vstart[1] = 7.67079e-05;
  vstart[2] = 3.64361e-01;

  step[0] = 2.69495e-03;
  step[1] = 8.78132e-10;
  step[2] = 2.35894e-06;
  gMinuit->mnparm(0, "a0", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
  gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);

  arglist[0] = 500;
  arglist[1] = 1.;
  cout << "Starting fitting" << endl;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  gMinuit->GetParameter(0, val1, err1);
  gMinuit->GetParameter(1, val2, err2);
  gMinuit->GetParameter(2, val3, err3);

  gMinuit ->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  cout << "val1: " << val1 << ", err1: " << err1 << endl;
  cout << "val2: " << val2 << ", err2: " << err2 << endl;
  cout << "val3: " << val3 << ", err3: " << err3 << endl;
  cout << "chi2: " << amin << ", ndf: " << size(Q2_values) - 3 << ", chi2/ndf: " << amin/(size(Q2_values) - 3) << endl;


  
  return 0;
}
