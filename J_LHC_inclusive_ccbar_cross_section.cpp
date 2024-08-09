#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <thread>
using namespace std;

#include "direct_dipole_amp_reader.h"

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

static const double Q2 = 0;

const double normalization = 8/M_PI*alpha_em*N_c*e_f*e_f;

//35, 34, 33, 32, 30, 25, 20, 15, 10, 5, 4, 3, 2, 1, 0.5
// 17, 16, 15, 10, 5, 1

static double r_limit; // 34.64
static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;
const string dipole_amp_type = "bfkl";
const string nucleus_type = "Pb";
const string filename_end = "";

const int warmup_calls = 100000;
const int integration_calls = 1000000;
const int integration_iterations = 1;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double x) {
  return get_dipole_amplitude(table, r, b_min, phi, x);
}

double L_integrand(double r, double b_min, double phi, double z, double Q2, double W) {
  return r*b_min*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, b_min, phi, (Q2+4*m_f*m_f)/(W*W+Q2));
}

double T_integrand(double r, double b_min, double phi, double z, double Q2, double W) {
  return r*b_min*(m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, b_min, phi, (Q2+4*m_f*m_f)/(W*W+Q2));
}

struct parameters {double W;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], Q2, par->W);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], Q2, par->W);
}

struct thread_par_struct
{
  double W;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double &a2, double &a3) : W(a1), sigma(a2), sigma_error(a3) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  const int dim = 4;
  double res, err;

  double xl[dim] = {0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1};
  params.W = par.W;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function L_G = {&L_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
  if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
    if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at W=" << params.W << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "L, W=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {


  const int dim = 4;
  double res, err;

  double xl[dim] = {0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1};
  params.W = par.W;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at W=" << params.W << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "T, W=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
  load_dipole_amplitudes(table, filename);

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

  /*
  double W = 100;
  double sigma_value;
  double sigma_error;
  thread_par_struct par(W, sigma_value, sigma_error);
  integrate_for_T_sigma(par);
  cout << "sigma=" << sigma_value << ", sigma_error=" << sigma_error << endl;
  return 0;
  */

  const int W_steps = 50;
  const double W_start = 2e1;
  const double W_stop = 2e4;
  const double W_step = 1.0/(W_steps-1)*log10(W_stop/W_start);

  stringstream r_limit_stream;
  r_limit_stream << fixed << setprecision(0) << r_limit;
  string r_limit_string = r_limit_stream.str();

  //stringstream r_limit_filename_stream;
  //r_limit_filename_stream << fixed << setprecision(0) << r_limit;
  //TString r_limit_filename_string = r_limit_filename_stream.str() + "_" +to_string(r_limit)[r_limit_filename_stream.str().length()+1]+to_string(r_limit)[r_limit_filename_stream.str().length()+2];
  TString r_limit_filename_string = r_limit_string;

  stringstream b_limit_stream;
  b_limit_stream << fixed << setprecision(0) << b_min_limit;
  string b_limit_string = b_limit_stream.str();

  TString b_limit_filename_string = b_limit_string;



  TMultiGraph* L_graphs = new TMultiGraph();
  TString title;
  if (print_r_limit) {
    title = "Longitudinal "+dipole_amp_type+" "+nucleus_type+" inclusive cross section with r limit: "+r_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else if (print_b_min_limit) {
    title = "Longitudinal "+dipole_amp_type+" "+nucleus_type+" inclusive cross section with b limit: "+b_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else {
    title = "Longitudinal "+dipole_amp_type+" "+nucleus_type+" inclusive cross section;W (GeV);cross section (mb)";
  }
  L_graphs->SetTitle(title);

  TString outfile_name;
  if (print_r_limit) {
    outfile_name = "data/J_LHC_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".txt";
  } else if (print_b_min_limit) {
    outfile_name = "data/J_LHC_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".txt";
  } else {
    outfile_name = "data/J_LHC_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+filename_end+".txt";
  }

  ofstream L_output_file(outfile_name);
  L_output_file << "W;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting L integration" << endl;
  double L_W_values[W_steps], L_sigma_values[W_steps], L_W_errors[W_steps], L_sigma_errors[W_steps];
  thread L_threads[W_steps];

  for (int i=0; i<W_steps; i++) {
    double W = pow(10, log10(W_start) + i*W_step);
    L_W_values[i] = W;
    L_W_errors[i] = 0;
    thread_par_struct par(W, L_sigma_values[i], L_sigma_errors[i]);
    L_threads[i] = thread(integrate_for_L_sigma, par);
  }

  for (int j=0; j<W_steps; j++) {
    L_threads[j].join();
  }
  TGraphErrors* subgraph = new TGraphErrors(W_steps, L_W_values, L_sigma_values, L_W_errors, L_sigma_errors);
  L_graphs->Add(subgraph);

  for (int i=0; i<W_steps; i++) {
    ostringstream W;
    W << L_W_values[i];
    ostringstream sigma;
    sigma << L_sigma_values[i];
    ostringstream sigma_err;
    sigma_err << L_sigma_errors[i];
    string line = W.str() + ";" + sigma.str() + ";" + sigma_err.str();
    L_output_file << line << endl;
  }
  
  L_output_file.close();

  TCanvas* L_sigma_canvas = new TCanvas("J_LHC_L_inclusive_sigma_canvas", "", 1100, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  gPad->SetLogy();

  //L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);


  if (print_r_limit) {
    outfile_name = "figures/J_LHC_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".pdf";
  } else if (print_b_min_limit) {
    outfile_name = "figures/J_LHC_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".pdf";
  } else {
    outfile_name = "figures/J_LHC_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+filename_end+".pdf";
  }
  L_sigma_canvas->Print(outfile_name);
 
 



  TMultiGraph* T_graphs = new TMultiGraph();

  if (print_r_limit) {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+" cross section with r limit: "+r_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else if (print_b_min_limit) {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+" cross section with b limit: "+b_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+" cross section;W (GeV);cross section (mb)";
  }
  T_graphs->SetTitle(title);

  if (print_r_limit) {
    outfile_name = "data/J_LHC_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".txt";
  } else if (print_b_min_limit) {
    outfile_name = "data/J_LHC_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".txt";
  } else {
    outfile_name = "data/J_LHC_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+filename_end+".txt";
  }
  ofstream T_output_file(outfile_name);
  T_output_file << "W;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting T integration" << endl;
  
  double T_W_values[W_steps], T_sigma_values[W_steps], T_W_errors[W_steps], T_sigma_errors[W_steps];
  thread T_threads[W_steps];

  for (int i=0; i<W_steps; i++) {
    double W = pow(10, log10(W_start) + i*W_step);
    T_W_values[i] = W;
    T_W_errors[i] = 0;
    thread_par_struct par(W, T_sigma_values[i], T_sigma_errors[i]);
    T_threads[i] = thread(integrate_for_T_sigma, par);
    //this_thread::sleep_for(30s);
  }

  for (int j=0; j<W_steps; j++) {
    T_threads[j].join();
  }

  subgraph = new TGraphErrors(W_steps, T_W_values, T_sigma_values, T_W_errors, T_sigma_errors);
  T_graphs->Add(subgraph);

  for (int i=0; i<W_steps; i++) {
    ostringstream W;
    W << T_W_values[i];
    ostringstream sigma;
    sigma << T_sigma_values[i];
    ostringstream sigma_err;
    sigma_err << T_sigma_errors[i];
    string line = W.str() + ";" + sigma.str() + ";" + sigma_err.str();
    T_output_file << line << endl;
  }
  
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("J_LHC_T_inclusive_sigma_canvas", "", 1100, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  gPad->SetLogy();

  //T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  if (print_r_limit) {
    outfile_name = "figures/J_LHC_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".pdf";
  } else if (print_b_min_limit) {
    outfile_name = "figures/J_LHC_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".pdf";
  } else {
    outfile_name = "figures/J_LHC_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+filename_end+".pdf";
  }
  T_sigma_canvas->Print(outfile_name);

  return 0;
}
