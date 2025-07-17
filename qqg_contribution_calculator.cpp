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
#include <ostream>
#include <sstream>
using namespace std;

#include "direct_dipole_amp_reader.h"

const double alpha_em = 1.0/137;
const int N_c = 3;
const double alpha_s = 0.25;

static double e_f = 2.0/3;
static double m_f = 1.27;

static double r_limit = 100; // 100
//static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

const int warmup_calls = 10000;
const int integration_calls = 100000;//20 000 000
const int integration_iterations = 1;

const string dipole_amp_type = "bk";
const string nucleus_type = "p";
string filename_end = "_all";
bool diffraction_dipamp = true;

const int i_start = 0; // number of data points to skip
const int data_inclusion_count = 226;

const int debug_precision = 10;
const double max_theta_root_excess = 1e-6;


const double sigma_0 = 70.26; //GeV^-2
const double C_f = 1.0;
const double Cyrille_x_0 = 1.632e-5;
const double lambda = 0.2197;
const double N_0 = 0.7;
const double gamma_c = 0.7376;
const double kappa = 9.9;
const double const_alpha = N_0*gamma_c/(8*(1-N_0));
const double const_beta = 1/2*exp(-(1-N_0)/(N_0*gamma_c)*log(1-N_0));


//static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
//static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

//static InterpMultilinear<4, double>* interpolator;

void read_data_file(string filename, vector<double> &Q2_values, vector<double> &beta_values, vector<double> &x_values, vector<double> &x_pom_F2_values, vector<double> &delta_stat_values, vector<double> &delta_sys_values) {
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
    beta_values.push_back(stod(value));
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
    x_pom_F2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    delta_stat_values.push_back(stod(value));
    i++;

    value = "";
    while(i < line.size()) {
      value += line[i];
      i++;
    }
    delta_sys_values.push_back(stod(value));
    i++;

  }
  cout << "Finished reading file" << endl;
}

double Q_s(double x) {
    return pow(Cyrille_x_0/x, lambda/2);
}

double no_b_dipamp(double rQ_s, double x) {

if (rQ_s <= 2) {
    return N_0*pow(rQ_s/2, 2*gamma_c)*exp(-1*(2*pow(log(rQ_s/2), 2))/(kappa*lambda*log(1/x)));
} else {
    return 1-exp(-4*const_alpha*pow(log(const_beta*rQ_s), 2));
}

}

struct Ig_parameters {
  double beta;
  double xpom;
  double Q2;
  double k2;
  double z;
};

double Ig_integrand(double r, void * parameters) {
  struct Ig_parameters * params = (struct Ig_parameters *)parameters;
  double beta = params->beta;
  double xpom = params->xpom;
  //double Q2 = params->Q2;
  double k2 = params->k2;
  double z = params->z;

  double x = beta*xpom;
  double a = sqrt(1-z);
  double b = sqrt(z);
  double c = Q_s(x)/sqrt(k2);
  return r*gsl_sf_bessel_Jn(2, a*r)*gsl_sf_bessel_Kn(2, b*r)*no_b_dipamp(c*r, x);
}

double Ig(double beta, double xpom, double Q2, double k2, double z) {

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  struct Ig_parameters parameters = {beta, xpom, Q2, k2, z};

  double result, error;

  gsl_function F;
  F.function = &Ig_integrand;
  F.params = &parameters;

  gsl_integration_qagiu(&F, 0, 0, 0.01, 100, w, &result, &error);

  gsl_integration_workspace_free(w);

  //cout << "Result: " << result << ", error: " << error << endl;

  return result;
}

double Fqqg_LLQ2_integrand(double beta, double xpom, double Q2, double k2, double z) {
  double normalization = sigma_0*alpha_s*C_f*N_c*beta*e_f*e_f/(32*gsl_pow_4(M_PI));
  return normalization*log(Q2/k2)*(gsl_pow_2(1-beta/z) + gsl_pow_2(beta/z))*gsl_pow_2(Ig(beta, xpom, Q2, k2, z));
}

struct qqg_LLQ2_parameters {
  double beta;
  double xpom;
  double Q2;
  double sigma;
  double sigma_error;
  double sigma_fit;
};

double integration_function_qqg_LLQ2(double *k, size_t dim, void * params) {

    double k2 = k[0];
    double z = k[1];
    struct qqg_LLQ2_parameters *par = (struct qqg_LLQ2_parameters *)params;

    return Fqqg_LLQ2_integrand(par->beta, par->xpom, par->Q2, k2, z);
}

double xpomFqqg_LLQ2(double beta, double xpom, double Q2, double &result, double &error, double &fit) {

  const int dim = 2;
  double res, err;

  double xl[dim] = {0, beta};
  double xu[dim] = {Q2, 1};

  struct qqg_LLQ2_parameters params = {1, 1, 1};
  params.beta = beta;
  params.xpom = xpom;
  params.Q2 = Q2;

  const gsl_rng_type *qqg_LLQ2_rng;
  gsl_rng *rng;

  gsl_monte_function qqg_LLQ2_rng_G = {&integration_function_qqg_LLQ2, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  qqg_LLQ2_rng = gsl_rng_default;
  rng = gsl_rng_alloc(qqg_LLQ2_rng);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&qqg_LLQ2_rng_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "qqg_LLQ2_warmup_error: " << status << endl; throw (status);}
  status = gsl_monte_vegas_integrate(&qqg_LLQ2_rng_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "qqg_LLQ2_integration_error: " << status << endl; throw (status);}

  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at xpom=" << params.xpom << endl;
  }
  result = res;
  error = err;
  fit = gsl_monte_vegas_chisq(T_s);
  cout << "Q²=" << params.Q2 << ", xpom=" << params.xpom << ", beta=" << params.beta << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);

  return 0;
}

double Fqqg_LLbeta_integrand(double beta, double xpom, double Q2, double r, double z, double Rx, double Ry) {
  double normalization = C_f*alpha_s*Q2*sigma_0/(8*gsl_pow_3(M_PI)*alpha_em);

  double epsilon = sqrt(z*(1-z)*Q2 + m_f*m_f);
  double Phi_T = alpha_em*N_c/(2*M_PI*M_PI)*e_f*e_f*((z*z+gsl_pow_2(1-z))*epsilon*epsilon*gsl_pow_2(gsl_sf_bessel_K1(epsilon*r)) + m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon*r)));

  double x = beta*xpom;
  double Qs = Q_s(x);
  double R = sqrt(Rx*Rx + Ry*Ry);
  double rminusR = sqrt(gsl_pow_2(r-Rx) + gsl_pow_2(0-Ry));
  double A = r*r/(R*R*gsl_pow_2(rminusR))*gsl_pow_2(no_b_dipamp(R*Qs, xpom) + no_b_dipamp(rminusR*Qs, xpom) - no_b_dipamp(r*Qs, xpom) - no_b_dipamp(R*Qs, xpom)*no_b_dipamp(rminusR*Qs, xpom));

  return normalization*r*Phi_T*A;
}

struct qqg_LLbeta_parameters {
  double beta;
  double xpom;
  double Q2;
  double sigma;
  double sigma_error;
  double sigma_fit;
};

double integration_function_qqg_LLbeta(double *k, size_t dim, void * params) {
  double beta = k[0];
  double r = k[1];
  double z = k[2];
  double Rx = k[3];
  double Ry = k[4];
  struct qqg_LLbeta_parameters *par = (struct qqg_LLbeta_parameters *)params;

  return Fqqg_LLbeta_integrand(par->xpom, par->Q2, beta, r, z, Rx, Ry);
}

double xpomFqqg_LLbeta(double beta, double xpom, double Q2, double &result, double &error, double &fit) {

  const int dim = 4;
  double res, err;

  double xl[dim] = {0, 0, -r_limit, -r_limit};
  double xu[dim] = {r_limit, 1, r_limit, r_limit};

  struct qqg_LLbeta_parameters params = {1, 1, 1};
  params.beta = beta;
  params.xpom = xpom;
  params.Q2 = Q2;

  const gsl_rng_type *qqg_LLbeta_rng;
  gsl_rng *rng;

  gsl_monte_function qqg_LLbeta_rng_G = {&integration_function_qqg_LLbeta, dim, &params};

  gsl_rng_env_setup();
  int status = 0;

  qqg_LLbeta_rng = gsl_rng_default;
  rng = gsl_rng_alloc(qqg_LLbeta_rng);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&qqg_LLbeta_rng_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "qqg_LLQ2_warmup_error: " << status << endl; throw (status);}
  status = gsl_monte_vegas_integrate(&qqg_LLbeta_rng_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "qqg_LLQ2_integration_error: " << status << endl; throw (status);}

  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at xpom=" << params.xpom << endl;
  }
  result = res;
  error = err;
  fit = gsl_monte_vegas_chisq(T_s);
  cout << "Q²=" << params.Q2 << ", xpom=" << params.xpom << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);

  return 0;
}

int main() {

  //gsl_set_error_handler_off();

  double result, error, fit;
  xpomFqqg_LLQ2(0.5, 3e-5, 1000, result, error, fit);
  cout << "result: " << result << ", error: " << error << ", fit: " << fit << endl << endl;

  xpomFqqg_LLbeta(0.5, 3e-5, 1000, result, error, fit);
  cout << "result: " << result << ", error: " << error << ", fit: " << fit << endl;

  return 0;

  vector<double> Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values;

  read_data_file("data/differential_HERA_data.dat", Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values);

  thread L_integration_threads[data_inclusion_count], T_integration_threads[data_inclusion_count];
  double L_sigma[data_inclusion_count], L_error[data_inclusion_count], L_fit[data_inclusion_count];
  double T_sigma[data_inclusion_count], T_error[data_inclusion_count], T_fit[data_inclusion_count];

  static auto t1 = chrono::high_resolution_clock::now();

  for (int i=0; i<data_inclusion_count; i++) {

    //thread_par_struct L_par(Q2_values[i+i_start], x_values[i+i_start], beta_values[i+i_start], L_sigma[i], L_error[i], L_fit[i]);
    //L_integration_threads[i] = thread(integrate_for_L_sigma, L_par);

    //thread_par_struct T_par(Q2_values[i+i_start], x_values[i+i_start], beta_values[i+i_start], T_sigma[i], T_error[i], T_fit[i]);
    //T_integration_threads[i] = thread(integrate_for_T_qqg, T_par);

  }

  for (int i=0; i<data_inclusion_count; i++) {
    //L_integration_threads[i].join();
    T_integration_threads[i].join();
  }

  static auto t2 = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
  cout << "Calculation finished in " << duration.count() << " seconds" << endl;

  ofstream L_output_file("output/differential_diffractive_L"+filename_end+".txt");
  L_output_file << "Q2 (GeV);beta;x;sigma (mb);sigma error (mb);fit" << endl;

  for (int i=0; i<data_inclusion_count; i++) {
    ostringstream Q2;
    Q2 << Q2_values[i+i_start];
    ostringstream beta;
    beta << beta_values[i+i_start];
    ostringstream x;
    x << x_values[i+i_start];
    ostringstream sigma;
    sigma << L_sigma[i];
    ostringstream error;
    error << L_error[i];
    ostringstream fit;
    fit << L_fit[i];
    string line = Q2.str() + ";" + beta.str() + ";" + x.str() + ";" + sigma.str() + ";" + error.str() + ";" + fit.str();
    L_output_file << line << endl;
  }
  L_output_file.close();

  ofstream T_output_file("output/differential_diffractive_T"+filename_end+".txt");
  T_output_file << "Q2 (GeV);beta;x;sigma (mb);sigma error (mb);fit" << endl;

  for (int i=0; i<data_inclusion_count; i++) {
    ostringstream Q2;
    Q2 << Q2_values[i+i_start];
    ostringstream beta;
    beta << beta_values[i+i_start];
    ostringstream x;
    x << x_values[i+i_start];
    ostringstream sigma;
    sigma << T_sigma[i];
    ostringstream error;
    error << T_error[i];
    ostringstream fit;
    fit << T_fit[i];
    string line = Q2.str() + ";" + beta.str() + ";" + x.str() + ";" + sigma.str() + ";" + error.str() + ";" + fit.str();
    T_output_file << line << endl;
  }
  T_output_file.close();
}
