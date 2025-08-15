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


const bool is_charm = true;

static double e_f;
static double m_f;

//static double r_limit = 100; // 100
//static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

const int warmup_calls = 1000;
const int integration_calls = 10000;//100000000
const int integration_iterations = 1;

const string dipole_amp_type = "bk";
const string nucleus_type = "p";
string filename_end = "_all";
bool diffraction_dipamp = true;

const int i_start = 0; // number of data points to skip
const int data_inclusion_count = 226;

const int debug_precision = 10;
const double max_theta_root_excess = 1e-6;


static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

static InterpMultilinear<4, double>* interpolator;


const double sigma_0 = 70.26; //GeV^-2
const double C_f = 1.0;
const double Cyrille_x_0 = 1.632e-5;
const double lambda = 0.2197;
const double N_0 = 0.7;
const double gamma_c = 0.7376;
const double kappa = 9.9;
const double const_alpha = N_0*gamma_c/(8*(1-N_0));
const double const_beta = 1.0/2*exp(-(1-N_0)/(N_0*gamma_c)*log(1-N_0));


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

double dipole_amplitude(double r, double b_min, double phi, double xpom) {
  //phi potentially being too large is ignored intentionally

  array<double, 4> args = {log(r), log(b_min), phi, log(xpom)};
  return exp(interpolator->interp(args.begin()));

}

struct Ig_parameters {
  double beta;
  double xpom;
  double Q2;
  double k2;
  double z;
};

double Ig_integrand(double r, double b_min, double phi, double beta, double xpom, double k2, double z) {

  double x = beta*xpom;
  double a = sqrt(1-z);
  double b = sqrt(z);
  double c = Q_s(x)/sqrt(k2);
  double N = dipole_amplitude(r, b_min, phi, xpom);
  double value = 2*r*gsl_sf_bessel_Jn(2, a*r)*gsl_sf_bessel_Kn(2, b*r)*(2*N-N*N);

  //cout << gsl_sf_bessel_Jn(2, a*r) << endl;
  //cout << gsl_sf_bessel_Kn(2, b*r) << endl;
  //cout << no_b_dipamp(c*r, x) << endl;
  //cout << "k2=" << k2 << endl;
  //cout << "beta=" << beta << ", z=" << z << ", r=" << r << ", rQs=" << c*r << ", x=" << x <<  ", value=" << value << endl << endl;
  return value;
}

double Ig_integration_function(double *k, size_t dim, void * parameters) {
  struct Ig_parameters * params = (struct Ig_parameters *)parameters;
  double beta = params->beta;
  double xpom = params->xpom;
  double k2 = params->k2;
  double z = params->z;
  double r = k[0];
  double b_min = k[1];
  double phi = k[2];

  return Ig_integrand(r, b_min, phi, beta, xpom, k2, z);
}

double Ig(double beta, double xpom, double Q2, double k2, double z) {

  const int dim = 3;
  double res, err;

  double xl[dim] = {0, 0, 0};
  double xu[dim] = {3.464101615137755e+01, 1.732050807568877e+01, M_PI};

  struct Ig_parameters parameters = {beta, xpom, Q2, k2, z};

  const gsl_rng_type *Ig_rng;
  gsl_rng *rng;

  gsl_monte_function Ig_G = {&Ig_integration_function, dim, &parameters};

  gsl_rng_env_setup ();
  int status = 0;
  
  Ig_rng = gsl_rng_default;
  rng = gsl_rng_alloc(Ig_rng);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *memory = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&Ig_G, xl, xu, dim, 100, rng, memory, &res, &err);
  if (status != 0) {cout << "Ig integration error: " << status << endl; throw (status);}
  status = gsl_monte_vegas_integrate(&Ig_G, xl, xu, dim, 1000, rng, memory, &res, &err);
  if (status != 0) {cout << "Ig integration error: " << status << endl; throw (status);}

  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at xpom=" << xpom << endl;
  }
  gsl_monte_vegas_free(memory);
  return res;
  //error = err;
  //double fit = gsl_monte_vegas_chisq(T_s);
  //cout << "LLQ2 Q²=" << params.Q2 << ", xpom=" << params.xpom << ", beta=" << params.beta << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;
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
  cout << "LLQ2 Q²=" << params.Q2 << ", xpom=" << params.xpom << ", beta=" << params.beta << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);

  return 0;
}

///////////////////////////////////////////


struct Ig_parameters {
  double beta;
  double xpom;
  double Q2;
  double k2;
  double z;
};

double nulbeta_Ig_integrand(double r, double b_min, double phi, double beta, double xpom, double k2, double z) {

  double N = dipole_amplitude(r, b_min, phi, xpom);
  double value = 2*1/r*gsl_sf_bessel_Jn(2, sqrt(k2)*r)*(2*N - N*N);

  //cout << gsl_sf_bessel_Jn(2, a*r) << endl;
  //cout << gsl_sf_bessel_Kn(2, b*r) << endl;
  //cout << no_b_dipamp(c*r, x) << endl;
  //cout << "k2=" << k2 << endl;
  //cout << "beta=" << beta << ", z=" << z << ", r=" << r << ", rQs=" << c*r << ", x=" << x <<  ", value=" << value << endl << endl;
  return value;
}

double nulbeta_Ig_integration_function(double *k, size_t dim, void * parameters) {
  struct Ig_parameters * params = (struct Ig_parameters *)parameters;
  double beta = params->beta;
  double xpom = params->xpom;
  double k2 = params->k2;
  double z = params->z;
  double r = k[0];
  double b_min = k[1];
  double phi = k[2];

  return nulbeta_Ig_integrand(r, b_min, phi, beta, xpom, k2, z);
}

double nulbeta_Ig(double beta, double xpom, double Q2, double k2, double z) {

  const int dim = 3;
  double res, err;

  double xl[dim] = {0, 0, 0};
  double xu[dim] = {3.464101615137755e+01, 1.732050807568877e+01, M_PI};

  struct Ig_parameters parameters = {beta, xpom, Q2, k2, z};

  const gsl_rng_type *Ig_rng;
  gsl_rng *rng;

  gsl_monte_function Ig_G = {&nulbeta_Ig_integration_function, dim, &parameters};

  gsl_rng_env_setup ();
  int status = 0;
  
  Ig_rng = gsl_rng_default;
  rng = gsl_rng_alloc(Ig_rng);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *memory = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&Ig_G, xl, xu, dim, 100, rng, memory, &res, &err);
  if (status != 0) {cout << "Ig integration error: " << status << endl; throw (status);}
  status = gsl_monte_vegas_integrate(&Ig_G, xl, xu, dim, 1000, rng, memory, &res, &err);
  if (status != 0) {cout << "Ig integration error: " << status << endl; throw (status);}

  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at xpom=" << xpom << endl;
  }
  gsl_monte_vegas_free(memory);
  return res;
  //error = err;
  //double fit = gsl_monte_vegas_chisq(T_s);
  //cout << "LLQ2 Q²=" << params.Q2 << ", xpom=" << params.xpom << ", beta=" << params.beta << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;
}

double nulbeta_Fqqg_LLQ2_integrand(double beta, double xpom, double Q2, double k2, double z) {
  double normalization = sigma_0*alpha_s*C_f*N_c*e_f*e_f/(12*gsl_pow_4(M_PI));
  return normalization*log(Q2/k2)*gsl_pow_2(nulbeta_Ig(beta, xpom, Q2, k2, z));
}

double nulbeta_integration_function_qqg_LLQ2(double *k, size_t dim, void * params) {

    double k2 = k[0];
    double z = 0;
    struct qqg_LLQ2_parameters *par = (struct qqg_LLQ2_parameters *)params;

    return nulbeta_Fqqg_LLQ2_integrand(par->beta, par->xpom, par->Q2, k2, z);
}

double nulbeta_xpomFqqg_LLQ2(double beta, double xpom, double Q2, double &result, double &error, double &fit) {

  const int dim = 1;
  double res, err;

  double xl[dim] = {0};
  double xu[dim] = {Q2};

  struct qqg_LLQ2_parameters params = {1, 1, 1};
  params.beta = beta;
  params.xpom = xpom;
  params.Q2 = Q2;

  const gsl_rng_type *qqg_LLQ2_rng;
  gsl_rng *rng;

  gsl_monte_function qqg_LLQ2_rng_G = {&nulbeta_integration_function_qqg_LLQ2, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  qqg_LLQ2_rng = gsl_rng_default;
  rng = gsl_rng_alloc(qqg_LLQ2_rng);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&qqg_LLQ2_rng_G, xl, xu, dim, 10, rng, T_s, &res, &err);
  if (status != 0) {cout << "qqg_LLQ2_warmup_error: " << status << endl; throw (status);}
  status = gsl_monte_vegas_integrate(&qqg_LLQ2_rng_G, xl, xu, dim, 100, rng, T_s, &res, &err);
  if (status != 0) {cout << "qqg_LLQ2_integration_error: " << status << endl; throw (status);}

  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at xpom=" << params.xpom << endl;
  }
  result = res;
  error = err;
  fit = gsl_monte_vegas_chisq(T_s);
  cout << "LLQ2 Q²=" << params.Q2 << ", xpom=" << params.xpom << ", beta=" << params.beta << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);

  return 0;
}

//////////////////////////////////////////

double polar_Fqqg_LLbeta_integrand(double beta, double xpom, double Q2, double r, double z, double R, double theta, double b_min, double phi) {
  double normalization = C_f*alpha_s*Q2*sigma_0/(8*gsl_pow_3(M_PI)*alpha_em);

  double epsilon = sqrt(z*(1-z)*Q2 + m_f*m_f);
  double Phi_T = alpha_em*N_c/(2*M_PI*M_PI)*e_f*e_f*((z*z+gsl_pow_2(1-z))*epsilon*epsilon*gsl_pow_2(gsl_sf_bessel_K1(epsilon*r)) + m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon*r)));

  //double x = beta*xpom;
  double Qs = Q_s(xpom);
  double rminusR = sqrt(r*r - 2*r*R*cos(theta) + R*R);
  double A = 2*2*R*r*r/(R*R*gsl_pow_2(rminusR))*gsl_pow_2(dipole_amplitude(R*Qs, b_min, phi, xpom) + dipole_amplitude(rminusR*Qs, b_min, phi, xpom) - dipole_amplitude(r*Qs, b_min, phi, xpom) - dipole_amplitude(R*Qs, b_min, phi, xpom)*dipole_amplitude(rminusR*Qs, b_min, phi, xpom));
  if (false) {
    cout << "beta=" << beta << ", xpom=" << xpom << ", Q2=" << Q2 << endl;
    cout << "Q_s: " << Qs << ", r: " << r << ", z: " << z << ", R: " << R << ", theta: " << theta << endl;
    cout << normalization*r*Phi_T*A << ", " << A << endl;;
    //cout << no_b_dipamp(R*Qs, xpom) << ", " << no_b_dipamp(rminusR*Qs, xpom) << ", -" << no_b_dipamp(r*Qs, xpom) << ", -" << no_b_dipamp(R*Qs, xpom)*no_b_dipamp(rminusR*Qs, xpom) << endl << endl;
  }
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
  double r = k[0];
  double z = k[1];
  double R = k[2];
  double theta = k[3];
  double b_min = k[4];
  double phi = k[5];
  struct qqg_LLbeta_parameters *par = (struct qqg_LLbeta_parameters *)params;

  return polar_Fqqg_LLbeta_integrand(par->beta, par->xpom, par->Q2, r, z, R, theta, b_min, phi);
}

double xpomFqqg_LLbeta(double beta, double xpom, double Q2, double &result, double &error, double &fit) {

  const int dim = 6;
  double res, err;

  //const double range = 10000;
  double xl[dim] = {0, 0, 0, 0, 0, 0};
  double xu[dim] = {100, 1, 100, M_PI, 100, M_PI};

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
  cout << "LLbeta: Q²=" << params.Q2 << ", xpom=" << params.xpom << ", beta=" << params.beta << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);

  return 0;
}

int calc_total_xpomF_Tqqg_contribution(double beta, double xpom, double Q2, double &result, double &error, double &fit) {
  double fits[3];

  double LLQ2_result, LLQ2_error, LLQ2_fit;
  xpomFqqg_LLQ2(beta, xpom, Q2, LLQ2_result, LLQ2_error, LLQ2_fit);

  double LLbeta_result, LLbeta_error, LLbeta_fit;
  xpomFqqg_LLbeta(beta, xpom, Q2, LLbeta_result, LLbeta_error, LLbeta_fit);

  double LLQ2nul_result, LLQ2nul_error, LLQ2nul_fit;
  nulbeta_xpomFqqg_LLQ2(beta, xpom, Q2, LLQ2nul_result, LLQ2nul_error, LLQ2nul_fit);

  result  = LLQ2_result*LLbeta_result/LLQ2nul_result;
  error = sqrt(gsl_pow_2(LLbeta_result/LLQ2nul_result*LLQ2_error) + gsl_pow_2(LLQ2_result/LLQ2nul_result*LLbeta_error) + gsl_pow_2(LLQ2_result*LLbeta_result/(LLQ2nul_result*LLQ2nul_result)*LLQ2nul_error));
  fits[0] = LLQ2_fit;
  fits[1] = LLbeta_fit;
  fits[2] = LLQ2nul_fit;
  if (LLQ2_fit > LLbeta_fit) {
    if (LLQ2_fit > LLQ2nul_fit) {
      fit = LLQ2_fit;
    } else {
      fit = LLQ2nul_fit;
    }
  } else {
    if (LLbeta_fit > LLQ2nul_fit) {
      fit = LLbeta_fit;
    } else {
      fit = LLQ2nul_fit;
    }
  }
  return 0;
}

struct thread_par_struct
{
  double Q2;
  double x;
  double beta;
  double &result;
  double &error;
  double &fit;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5, double &a6) : Q2(a1), x(a2), beta(a3), result(a4), error(a5), fit(a6) {}
};

int low_beta_thread_func(thread_par_struct par) {
  calc_total_xpomF_Tqqg_contribution(par.beta, par.x, par.Q2, par.result, par.error, par.fit);
  return 0;
}

int main() {
  if (is_charm) {
    e_f = 2.0/3;
    m_f = 1.27;
    filename_end += "_charm";
  } else {
    e_f = sqrt(2.0/3*2.0/3+2*1.0/3*1.0/3);
    m_f = 0;
  }

  gsl_set_error_handler_off();

  double beta = 0.04;
  double xpom = 0.00012/0.04;
  double Q2 = 4.5;
  double r = 1e-3;
  cout << r << ", " << polar_Fqqg_LLbeta_integrand(beta, xpom, Q2, r, 0.5, r, 0.01) << endl;;
  return 0;
  /*
  double aresult, aerror, afit;
  calc_total_xpomF_Tqqg_contribution(0.04, 0.00012/0.04, 4.5, aresult, aerror, afit);
  cout << aresult << endl;
  cout << aerror << endl;
  cout << afit << endl;
  return 0;
  */
  /*
  for(int i=0; i<1000; i++) {
    double beta = 0.04;
    double xpom = 0.00012/0.04;
    double Q2 = 4.5;
    double r = i/1000.0*10+1.0/1000;
    cout << r << ", " << polar_Fqqg_LLbeta_integrand(beta, xpom, Q2, r, 0.5, r, 0.01) << endl;;
  }
  return 0;
  */

  //gsl_set_error_handler_off();
  /*
  double Result, Error, Fit;

  xpomFqqg_LLbeta(0.04, 0.00012/0.04, 4.5, Result, Error, Fit);
  cout << "result: " << Result << ", error: " << Error << ", fit: " << Fit << endl;
  return 0;
  */
  vector<double> Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values;

  read_data_file("data/differential_HERA_data.dat", Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values);

  //string filename = "data/dipole_amplitude_with_IP_dependence_bk_p_diffraction.csv";
  //load_p_dipole_amplitudes(table, filename);

  thread low_beta_integration_threads[data_inclusion_count];
  double result[data_inclusion_count], error[data_inclusion_count], fit[data_inclusion_count];

  static auto t1 = chrono::high_resolution_clock::now();

  for (int i=0; i<data_inclusion_count; i++) {

    thread_par_struct low_beta_par(Q2_values[i+i_start], x_values[i+i_start]/beta_values[i+i_start], beta_values[i+i_start], result[i], error[i], fit[i]);
    low_beta_integration_threads[i] = thread(low_beta_thread_func, low_beta_par);

  }

  for (int i=0; i<data_inclusion_count; i++) {
    low_beta_integration_threads[i].join();
  }

  static auto t2 = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
  cout << "Calculation finished in " << duration.count() << " seconds" << endl;

  ofstream L_output_file("output/low_beta_corrections"+filename_end+".txt");
  L_output_file << "Q2 (GeV);beta;x;sigma (mb);sigma error (mb);fit" << endl;

  for (int i=0; i<data_inclusion_count; i++) {
    ostringstream Q2;
    Q2 << Q2_values[i+i_start];
    ostringstream beta;
    beta << beta_values[i+i_start];
    ostringstream x;
    x << x_values[i+i_start];
    ostringstream result_stream;
    result_stream << result[i];
    ostringstream error_stream;
    error_stream << error[i];
    ostringstream fit_stream;
    fit_stream << fit[i];
    string line = Q2.str() + ";" + beta.str() + ";" + x.str() + ";" + result_stream.str() + ";" + error_stream.str() + ";" + fit_stream.str();
    L_output_file << line << endl;
  }
  L_output_file.close();
}

