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
const double e_f = sqrt(2.0/3*2.0/3+1.0/3*1.0/3+1.0/3*1.0/3); //2.0/3, sqrt(2.0/3*2.0/3+1.0/3*1.0/3+1.0/3*1.0/3)
const double m_f = 0; //GeV 1.27

const double normalization = 1.0/gsl_pow_3(2*M_PI)*alpha_em*N_c*e_f*e_f;

static double r_limit; // 34.64
static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

const int warmup_calls = 100000;
const int integration_calls = 1000000;//20 000 000
const int integration_iterations = 1;

const string dipole_amp_type = "bfkl";
const string nucleus_type = "p";
const string filename_end = "";

const int i_start = 0; // number of data points to skip
const int data_inclusion_count = 20;

const string filename_end = "_direct_1mil_0-19_xpom";//

const int debug_precision = 10;
const double max_theta_root_excess = 1e-6;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

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

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double x) {
  if (nucleus_type == "p") {
    return get_p_dipole_amplitude(p_table, r, b_min, phi, x, false);
  } else if (nucleus_type == "Pb") {
    return get_Pb_dipole_amplitude(Pb_table, r, b_min, phi, x, false);
  } else {
    throw 1;
  }
}

double calc_phi(double r1, double r2, double b1, double b2, double z) {
  double phi;
  double raw_phi = atan2(-r2, -r1) - atan2((b2+(1-z)*r2), (b1+(1-z)*r1));
  if (M_PI < raw_phi) {
    phi = 2*M_PI - raw_phi;
  } else if (-M_PI <= raw_phi && raw_phi < 0) {
    phi = -raw_phi;
  } else if (raw_phi < -M_PI) {
    phi = 2*M_PI + raw_phi;
  } else {
    phi = raw_phi;
  }
  return phi;
}

double L_integrand(double r1, double r2, double b1, double b2, double r1bar, double r2bar, double z, double Q2, double x, double beta) {
  if (z*(1-z)*Q2*(1/beta-1)-m_f*m_f < 0) {
    cout << "setting zero, z=" << z << endl;
    return 0;
  }
  
  double r = sqrt(r1*r1 + r2*r2);
  double rbar = sqrt(r1bar*r1bar + r2bar*r2bar);
  double bmin = sqrt(gsl_pow_2(b1+(1-z)*r1) + gsl_pow_2(b2+(1-z)*r2));
  double bminbar = sqrt(gsl_pow_2(b1+(1-z)*r1bar) + gsl_pow_2(b2+(1-z)*r2bar));
  double phi = calc_phi(r1, r2, b1, b2, z);
  double phibar = calc_phi(r1bar, r2bar, b1, b2, z);

  return gsl_sf_bessel_J0(sqrt(z*(1-z)*Q2*(1/beta-1)-m_f*m_f)*sqrt(gsl_pow_2(r1-r1bar)+gsl_pow_2(r2-r2bar)))
  *z*(1-z)
  *4*Q2*z*z*gsl_pow_2(1-z)*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*rbar)
  *dipole_amplitude(r, bmin, phi, x/beta)*dipole_amplitude(rbar, bminbar, phibar, x/beta);

}

double T_integrand(double r1, double r2, double b1, double b2, double r1bar, double r2bar, double z, double Q2, double x, double beta) {
  if (z*(1-z)*Q2*(1/beta-1)-m_f*m_f < 0) {
    cout << "setting zero, z=" << z << endl;
    return 0;
  }
  
  double r = sqrt(r1*r1 + r2*r2);
  double rbar = sqrt(r1bar*r1bar + r2bar*r2bar);
  double bmin = sqrt(gsl_pow_2(b1+(1-z)*r1) + gsl_pow_2(b2+(1-z)*r2));
  double bminbar = sqrt(gsl_pow_2(b1+(1-z)*r1bar) + gsl_pow_2(b2+(1-z)*r2bar));
  double phi = calc_phi(r1, r2, b1, b2, z);
  double phibar = calc_phi(r1bar, r2bar, b1, b2, z);

  return gsl_sf_bessel_J0(sqrt(z*(1-z)*Q2*(1/beta-1)-m_f*m_f)*sqrt(gsl_pow_2(r1-r1bar)+gsl_pow_2(r2-r2bar)))
  *z*(1-z)
  *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*rbar) + epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*(r1*r1bar+r2*r2bar)/(r*rbar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*rbar))
  *dipole_amplitude(r, bmin, phi, x/beta)*dipole_amplitude(rbar, bminbar, phibar, x/beta);
}

struct parameters {double Q2; double x; double beta;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], k[4], k[5], k[6], par->Q2, par->x, par->beta);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], k[6], par->Q2, par->x, par->beta);
}

struct thread_par_struct
{
  double Q2;
  double x;
  double beta;
  double &sigma;
  double &sigma_error;
  double &sigma_fit;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5, double &a6) : Q2(a1), x(a2), beta(a3), sigma(a4), sigma_error(a5), sigma_fit(a6) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

    const int dim = 7;
    double res, err;

    double z_min = (1-sqrt(1-4*m_f*m_f/(par.Q2*(1/par.beta-1))))/2;
    double z_max = (1+sqrt(1-4*m_f*m_f/(par.Q2*(1/par.beta-1))))/2;

    double xl[dim] = {-r_limit, -r_limit, -b_min_limit, -b_min_limit, -r_limit, -r_limit, z_min};
    double xu[dim] = {r_limit, r_limit, b_min_limit, b_min_limit, r_limit, r_limit, z_max};

    struct parameters params = {1, 1, 1};
    params.Q2 = par.Q2;
    params.x = par.x;
    params.beta = par.beta;
    double &sigma = par.sigma;
    double &sigma_error = par.sigma_error;
    double &sigma_fit = par.sigma_fit;

    const gsl_rng_type *T;
    gsl_rng *rng;

    gsl_monte_function L_G = {&L_g, dim, &params};

    gsl_rng_env_setup ();
    int status = 0;

    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, 1);

    gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);
    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
    if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
    for (int i=0; i<integration_iterations; i++) {
        static auto t1 = chrono::high_resolution_clock::now();
        status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
        if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
        static auto t2 = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
        //cout << "L iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(L_s) << ", duration: " << duration.count() << endl;
    }
    if (gsl_isnan(res)) {
        res = 0;
        cout << "nan found at x=" << params.x << endl;
    }
    sigma = res;
    sigma_error = err;
    sigma_fit = gsl_monte_vegas_chisq(L_s);
    cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

    gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 7;
  double res, err;

  double z_min = (1-sqrt(1-4*m_f*m_f/(par.Q2*(1/par.beta-1))))/2;
  double z_max = (1+sqrt(1-4*m_f*m_f/(par.Q2*(1/par.beta-1))))/2;

  double xl[dim] = {-r_limit, -r_limit, -b_min_limit, -b_min_limit, -r_limit, -r_limit, z_min};
  double xu[dim] = {r_limit, r_limit, b_min_limit, b_min_limit, r_limit, r_limit, z_max};

  struct parameters params = {1, 1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
  params.beta = par.beta;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;
  double &sigma_fit = par.sigma_fit;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    static auto t1 = chrono::high_resolution_clock::now();
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    //cout << "T iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(T_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  sigma_fit = gsl_monte_vegas_chisq(T_s);
  cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", beta=" << params.beta << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  vector<double> Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values;

  read_data_file("data/differential_HERA_data.dat", Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values);

  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
  if (nucleus_type == "p") {
    load_p_dipole_amplitudes(p_table, filename);
  } else if (nucleus_type == "Pb") {
    load_Pb_dipole_amplitudes(Pb_table, filename);
  } else {
    throw 1;
  }

  if (nucleus_type == "Pb") {
    r_limit = 657; // 34.64, 657
    b_min_limit = 328; // 17.32, 328
  } else if (nucleus_type == "p") {
    r_limit = 34.64;
    b_min_limit = 17.32;
  } else {
    cout << "invalid nucleus type" << endl;
    throw 1;
  }

  thread L_integration_threads[data_inclusion_count], T_integration_threads[data_inclusion_count];
  double L_sigma[data_inclusion_count], L_error[data_inclusion_count], L_fit[data_inclusion_count];
  double T_sigma[data_inclusion_count], T_error[data_inclusion_count], T_fit[data_inclusion_count];

  static auto t1 = chrono::high_resolution_clock::now();

  for (int i=0; i<data_inclusion_count; i++) {

    thread_par_struct L_par(Q2_values[i+i_start], x_values[i+i_start], beta_values[i+i_start], L_sigma[i], L_error[i], L_fit[i]);
    L_integration_threads[i] = thread(integrate_for_L_sigma, L_par);

    thread_par_struct T_par(Q2_values[i+i_start], x_values[i+i_start], beta_values[i+i_start], T_sigma[i], T_error[i], T_fit[i]);
    T_integration_threads[i] = thread(integrate_for_T_sigma, T_par);

  }

  for (int i=0; i<data_inclusion_count; i++) {
    L_integration_threads[i].join();
    T_integration_threads[i].join();
  }

  static auto t2 = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
  cout << "Calculation finished in " << duration.count() << " seconds" << endl;

  ofstream L_output_file("data/differential_diffractive_L"+filename_end+".txt");
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

  ofstream T_output_file("data/differential_diffractive_T"+filename_end+".txt");
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
