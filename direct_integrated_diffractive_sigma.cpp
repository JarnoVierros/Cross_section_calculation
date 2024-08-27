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
const double e_f = 2.0/3;
const double m_f = 1.27;

const double normalization = 2.0/gsl_pow_2(2*M_PI)*alpha_em*N_c*e_f*e_f;

static double r_limit; // 34.64
static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

static int warmup_calls;
static int integration_calls;
const int integration_iterations = 1;

//const string filename_end = "_20mil_85-225";//

const int debug_precision = 10;

const string dipole_amp_type = "bfkl";
const string nucleus_type = "p";
const string filename_end = "_special_1";

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double x, double Q2) {
  double shifted_x = x*(1+4*m_f*m_f/Q2);
  if (nucleus_type == "p") {
    return get_p_dipole_amplitude(p_table, r, b_min, phi, x, false);
  } else if (nucleus_type == "Pb") {
    return get_Pb_dipole_amplitude(Pb_table, r, b_min, phi, x, false);
  } else {
    throw 1;
  }
}

double calc_phi(double r1, double r2, double b1, double b2, double z) {
  
  double r = sqrt(r1*r1 + r2*r2);
  double b = sqrt(b1*b1 + b2*b2);
  double X = sqrt(gsl_pow_2(1-z)*r*r+2*(1-z)*(r1*b1+r2*b2)+b*b);
  double Y = sqrt(z*z*r*r-2*z*(b1*r1+b2*r2)+b*b);
  double phi;
  if (X <= Y) {
    double raw_phi = abs(atan2(r2, r1) - atan2(b2+(1-z)*r2, b1+(1-z)*r1) + M_PI);
    if (raw_phi > M_PI) {
      phi = 2*M_PI - raw_phi;
    } else {
      phi = raw_phi;
    }
  } else {
    double raw_phi = abs(atan2(-r2, -r1)-atan2(b2-z*r2, b1-z*r1)+M_PI);
    if (raw_phi > M_PI) {
      phi = 2*M_PI - raw_phi;
    } else {
      phi = raw_phi;
    }
  }
  phi = abs(phi);
  return phi;
}

double calc_bmin(double r1, double r2, double b1, double b2, double z) {
  double r = sqrt(r1*r1 + r2*r2);
  double b = sqrt(b1*b1 + b2*b2);
  double X = sqrt(gsl_pow_2(1-z)*r*r+2*(1-z)*(r1*b1+r2*b2)+b*b);
  double Y = sqrt(z*z*r*r-2*z*(b1*r1+b2*r2)+b*b);
  double bmin;
  if (X < Y) {
    bmin = sqrt(gsl_pow_2(b1 + (1-z)*r1) + gsl_pow_2(b2 + (1-z)*r2));
  } else {
    bmin = sqrt(gsl_pow_2(b1-z*r1) + gsl_pow_2(b2-z*r2));
  }
  return bmin;
}

double max_sus_integrand = 0;

double L_integrand(double r1, double r2, double b1, double b2, double z, double Q2, double x) {
  
  double r = sqrt(r1*r1 + r2*r2);
  double bmin = calc_bmin(r1, r2, b1, b2, z);
  double phi = calc_phi(r1, r2, b1, b2, z);

  return 4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))
  *gsl_pow_2(dipole_amplitude(r, bmin, phi, x, Q2));

}

double T_integrand(double r1, double r2, double b1, double b2, double z, double Q2, double x) {
  
  double r = sqrt(r1*r1 + r2*r2);
  double bmin = calc_bmin(r1, r2, b1, b2, z);
  double phi = calc_phi(r1, r2, b1, b2, z);

  return (m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) 
  + epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))
  *gsl_pow_2(dipole_amplitude(r, bmin, phi, x, Q2));
  
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], k[4], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], k[4], par->Q2, par->x);
}

struct thread_par_struct
{
  double Q2;
  double x;
  double &sigma;
  double &sigma_error;
  double &sigma_fit;
  int seed;
  thread_par_struct(double a1, double a2, double &a3, double &a4, double &a5, int a6) : Q2(a1), x(a2), sigma(a3), sigma_error(a4), sigma_fit(a5), seed(a6) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

    const int dim = 5;
    double res, err;

    double xl[dim] = {-r_limit, -r_limit, -b_min_limit, -b_min_limit, 0};
    double xu[dim] = {r_limit, r_limit, b_min_limit, b_min_limit, 1};

    struct parameters params = {1, 1};
    params.Q2 = par.Q2;
    params.x = par.x;
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
    gsl_rng_set(rng, par.seed);

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
    cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

    gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 5;
  double res, err;

  double xl[dim] = {-r_limit, -r_limit, -b_min_limit, -b_min_limit, 0};
  double xu[dim] = {r_limit, r_limit, b_min_limit, b_min_limit, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
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
  gsl_rng_set(rng, par.seed);

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
  cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

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
    warmup_calls = 100000;
    integration_calls = 5000000;
  } else if (nucleus_type == "p") {
    r_limit = 34.64;
    b_min_limit = 17.32;
    warmup_calls = 100000;
    integration_calls = 3000000;
  } else {
    cout << "invalid nucleus type" << endl;
    throw 1;
  }

  const int Q2_values[] = {1, 3, 5, 8, 10};

  const int x_steps = 30;
  const double x_start = 1e-5;
  const double x_stop = 0.01;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Diffractive longitudinal cross section;x;cross section (mb)");

  ofstream L_output_file("data/direct_diff_L_sigma_"+dipole_amp_type+"_"+nucleus_type+".txt");
  L_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting L integration" << endl;
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double L_x_values[x_steps], L_sigma_values[x_steps], L_x_errors[x_steps], L_sigma_errors[x_steps], L_sigma_fits[x_steps];
    thread L_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      L_x_values[i] = x;
      L_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, L_sigma_values[i], L_sigma_errors[i], L_sigma_fits[i], 0);
      L_threads[i] = thread(integrate_for_L_sigma, par);
    }

    for (int j=0; j<x_steps; j++) {
      L_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(x_steps, L_x_values, L_sigma_values, L_x_errors, L_sigma_errors);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);

    for (int i=0; i<x_steps; i++) {
      ostringstream x;
      x << L_x_values[i];
      ostringstream sigma;
      sigma << L_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << L_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + x.str() + ";" + sigma.str() + ";" + sigma_err.str()+ ";" + sigma_err.str();
      L_output_file << line << endl;
    }
  }
  L_output_file.close();

  TCanvas* L_sigma_canvas = new TCanvas("diff_L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  TString title = "figures/direct_diff_L_sigma_"+dipole_amp_type+"_"+nucleus_type+".txt";
  L_sigma_canvas->Print(title);
 

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Diffractive transverse cross section;x;cross section (mb)");

  ofstream T_output_file("data/direct_diff_T_sigma_"+dipole_amp_type+"_"+nucleus_type+".txt");
  T_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting T integration" << endl;
  
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps], T_sigma_fits[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i], 0);
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

    for (int i=0; i<x_steps; i++) {
      ostringstream x;
      x << T_x_values[i];
      ostringstream sigma;
      sigma << T_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << T_sigma_errors[i];
      string line = to_string(Q2_values[j]) + ";" + x.str() + ";" + sigma.str() + ";" + sigma_err.str();
      T_output_file << line << endl;
    }
  }
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  title = "figures/direct_diff_T_sigma_"+dipole_amp_type+"_"+nucleus_type+".pdf";
  T_sigma_canvas->Print(title);

  return 0;

}