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
static double e_f = 2.0/3;
static double m_f = 1.27; //GeV

const double normalization = 8/M_PI*alpha_em*N_c*e_f*e_f;

static double r_limit; // 34.64
static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;
const string dipole_amp_type = "bk";
const string nucleus_type = "p";
const string filename_end = "";

const int warmup_calls = 100000;
const int integration_calls = 1000000;
const int integration_iterations = 1;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

static InterpMultilinear<4, double>* interpolator;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double W, double Q2) {
  double shifted_x = (Q2+4*m_f*m_f)/(W*W+Q2);

  if (calc_max_phi(r, b_min) < phi) {
    return 0;
  } else {
    array<double, 4> args = {log(r), log(b_min), phi, log(shifted_x)};
    return exp(interpolator->interp(args.begin()));
  }
}

double L_integrand(double r, double b_min, double phi, double z, double Q2, double W) {
  return r*b_min*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, b_min, phi, W, Q2);
}

double T_integrand(double r, double b_min, double phi, double z, double Q2, double W) {
  return r*b_min*(m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, b_min, phi, W, Q2);
}

struct parameters {double Q2; double W;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], par->Q2, par->W);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], par->Q2, par->W);
}

struct thread_par_struct
{
  double Q2;
  double W;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double a2, double &a3, double &a4) : Q2(a1), W(a2), sigma(a3), sigma_error(a4) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  const int dim = 4;
  double res, err;

  double xl[dim] = {0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
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
    cout << "nan found at Q2=" << params.Q2 << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "L, Q²=" << params.Q2 << ", W=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {


  const int dim = 4;
  double res, err;

  double xl[dim] = {0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
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
  cout << "T, Q²=" << params.Q2 << ", W=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  const double W_values[] = {60, 90, 120}; //{60, 90, 120};

  const int Q2_steps = 20;
  const double Q2_start = 1e-1;
  const double Q2_stop = 1e2;
  const double Q2_step = 1.0/(Q2_steps-1)*log10(Q2_stop/Q2_start);

  stringstream r_limit_stream;
  r_limit_stream << fixed << setprecision(0) << r_limit;
  string r_limit_string = r_limit_stream.str();

  TString r_limit_filename_string = r_limit_string;

  stringstream b_limit_stream;
  b_limit_stream << fixed << setprecision(0) << b_min_limit;
  string b_limit_string = b_limit_stream.str();

  TString b_limit_filename_string = b_limit_string;

  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
  if (nucleus_type == "p") {
    load_p_dipole_amplitudes(p_table, filename);
    create_p_interpolator(p_table, interpolator);
  } else if (nucleus_type == "Pb") {
    load_Pb_dipole_amplitudes(Pb_table, filename);
    create_Pb_interpolator(Pb_table, interpolator);
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

  TMultiGraph* L_graphs = new TMultiGraph();
  TString title;
  if (print_r_limit) {
    title = "Longitudinal "+dipole_amp_type+" "+nucleus_type+" inclusive cross section with r limit: "+r_limit_string+" GeV^-1;x;cross section (mb)";
  } else if (print_b_min_limit) {
    title = "Longitudinal "+dipole_amp_type+" "+nucleus_type+" inclusive cross section with b limit: "+b_limit_string+" GeV^-1;x;cross section (mb)";
  } else {
    title = "Longitudinal "+dipole_amp_type+" "+nucleus_type+" inclusive cross section;Q^{2} (GeV^{2});cross section (GeV^{-2})";
  }
  L_graphs->SetTitle(title);

  TString outfile_name;
  if (print_r_limit) {
    outfile_name = "data/J_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".txt";
  } else if (print_b_min_limit) {
    outfile_name = "data/J_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".txt";
  } else {
    outfile_name = "data/J_L_inclusive_Q2_"+dipole_amp_type+"_"+nucleus_type+filename_end+".txt";
  }

  ofstream L_output_file(outfile_name);
  L_output_file << "W (GeV);Q2 (GeV^2);sigma (GeV^-2);sigma error (GeV^-2)" << endl;

  cout << "Starting L integration" << endl;
  for (long unsigned int j=0; j<size(W_values); j++) {
    double L_Q2_values[Q2_steps], L_sigma_values[Q2_steps], L_Q2_errors[Q2_steps], L_sigma_errors[Q2_steps];
    thread L_threads[Q2_steps];

    for (int i=0; i<Q2_steps; i++) {
      double Q2 = pow(10, log10(Q2_start) + i*Q2_step);
      L_Q2_values[i] = Q2;
      L_Q2_errors[i] = 0;
      thread_par_struct par(Q2, W_values[j], L_sigma_values[i], L_sigma_errors[i]);
      L_threads[i] = thread(integrate_for_L_sigma, par);
    }

    for (int j=0; j<Q2_steps; j++) {
      L_threads[j].join();
    }
    TGraphErrors* subgraph = new TGraphErrors(Q2_steps, L_Q2_values, L_sigma_values, L_Q2_errors, L_sigma_errors);

    stringstream stream;
    stream << fixed << setprecision(1) << W_values[j];
    TString subgraph_name = "W=" + stream.str() +" GeV";
    subgraph->SetTitle(subgraph_name);
    L_graphs->Add(subgraph);

    for (int i=0; i<Q2_steps; i++) {
      ostringstream Q2;
      Q2 << L_Q2_values[i];
      ostringstream sigma;
      sigma << L_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << L_sigma_errors[i];
      string line = to_string(W_values[j]) + ";" + Q2.str() + ";" + sigma.str() + ";" + sigma_err.str();
      L_output_file << line << endl;
    }
  }
  L_output_file.close();

  TCanvas* L_sigma_canvas = new TCanvas("J_L_inclusive_Q2_sigma_canvas", "", 1100, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);


  if (print_r_limit) {
    outfile_name = "figures/J_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".pdf";
  } else if (print_b_min_limit) {
    outfile_name = "figures/J_L_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".pdf";
  } else {
    outfile_name = "figures/J_L_inclusive_Q2_"+dipole_amp_type+"_"+nucleus_type+filename_end+".pdf";
  }
  L_sigma_canvas->Print(outfile_name);
 
 



  TMultiGraph* T_graphs = new TMultiGraph();

  if (print_r_limit) {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+" cross section with r limit: "+r_limit_string+" GeV^-1;x;cross section (mb)";
  } else if (print_b_min_limit) {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+" cross section with b limit: "+b_limit_string+" GeV^-1;x;cross section (mb)";
  } else {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+" cross section;x;cross section (mb)";
  }
  T_graphs->SetTitle(title);

  if (print_r_limit) {
    outfile_name = "data/J_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".txt";
  } else if (print_b_min_limit) {
    outfile_name = "data/J_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".txt";
  } else {
    outfile_name = "data/J_T_inclusive_Q2_"+dipole_amp_type+"_"+nucleus_type+filename_end+".txt";
  }
  ofstream T_output_file(outfile_name);
  T_output_file << "W (GeV);Q2 (GeV^2);sigma (GeV^-2);sigma error (GeV^-2)" << endl;

  cout << "Starting T integration" << endl;
  
  for (long unsigned int j=0; j<size(W_values); j++) {
    double T_Q2_values[Q2_steps], T_sigma_values[Q2_steps], T_Q2_errors[Q2_steps], T_sigma_errors[Q2_steps];
    thread T_threads[Q2_steps];

    for (int i=0; i<Q2_steps; i++) {
      double Q2 = pow(10, log10(Q2_start) + i*Q2_step);
      T_Q2_values[i] = Q2;
      T_Q2_errors[i] = 0;
      thread_par_struct par(Q2, W_values[j], T_sigma_values[i], T_sigma_errors[i]);
      T_threads[i] = thread(integrate_for_T_sigma, par);
      //this_thread::sleep_for(30s);
    }

    for (int j=0; j<Q2_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(Q2_steps, T_Q2_values, T_sigma_values, T_Q2_errors, T_sigma_errors);
    
    stringstream stream;
    stream << fixed << setprecision(1) << W_values[j];
    TString subgraph_name = "W=" + stream.str() + " GeV";
    subgraph->SetTitle(subgraph_name);
    T_graphs->Add(subgraph);

    for (int i=0; i<Q2_steps; i++) {
      ostringstream Q2;
      Q2 << T_Q2_values[i];
      ostringstream sigma;
      sigma << T_sigma_values[i];
      ostringstream sigma_err;
      sigma_err << T_sigma_errors[i];
      string line = to_string(W_values[j]) + ";" + Q2.str() + ";" + sigma.str() + ";" + sigma_err.str();
      T_output_file << line << endl;
    }
  }
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("J_T_inclusive_Q2_sigma_canvas", "", 1100, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  if (print_r_limit) {
    outfile_name = "figures/J_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".pdf";
  } else if (print_b_min_limit) {
    outfile_name = "figures/J_T_inclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".pdf";
  } else {
    outfile_name = "figures/J_T_inclusive_Q2_"+dipole_amp_type+"_"+nucleus_type+filename_end+".pdf";
  }
  T_sigma_canvas->Print(outfile_name);

  return 0;
}
