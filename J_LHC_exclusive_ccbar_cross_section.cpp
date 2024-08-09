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

const double normalization = 16.0/gsl_pow_2(2*M_PI)*alpha_em*N_c*e_f*e_f;



const bool print_r_limit = false;
const bool print_b_min_limit = false;

const int warmup_calls = 100000;
const int integration_calls = 20000000;//20 000 000
const int integration_iterations = 1;

const string dipole_amp_type = "bk";
const string nucleus_type = "Pb";
const string filename_end = "";

const double r_limit;
const double b_min_limit;

//const string filename_end = "_20mil_85-225";//

const int debug_precision = 10;
const double max_theta_root_excess = 1e-6;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

double calc_h(double r, double b_min, double phi, double z) {
  return b_min*b_min + (1-z)*2*b_min*r*cos(phi) + gsl_pow_2(1-z)*r*r;
}

double calc_b1(double r, double b_min, double phi, double z) {
  return b_min + (1-z)*r*cos(phi);
}

double calc_b2(double r, double b_min, double phi, double z) {
  return (1-z)*r*sin(phi);
}

double calc_j(double b2, double r_bar, double phi_bar, double z) {
  return (1-z)*2*b2*r_bar*sin(phi_bar);
}

double calc_A(double j, double h, double b1, double b2) {
  return sqrt(j*j + 4*h*(b1*b1-gsl_pow_2(j/(2*b2))));
}

bool theta_root_invalid(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar, double z) {
  double excess = abs(tan(theta_bar)*(b_min+(1-z)*r*cos(phi)-(1-z)*r_bar*cos(theta_bar+phi_bar)) + (1-z)*r_bar*sin(theta_bar+phi_bar) - (1-z)*r*sin(phi));
  if (excess > max_theta_root_excess) {
    return 1;
  } else {
    return 0;
  }
}

int calc_theta_bar(double return_values[4], double r, double b_min, double phi, double r_bar, double phi_bar, double z) {
  double b1 = calc_b1(r, b_min, phi, z);
  double b2 = calc_b2(r, b_min, phi, z);
  double h = calc_h(r, b_min, phi, z);
  if (r_bar*r_bar > (4*h*b1*b1)/(gsl_pow_2((1-z)*2*sin(phi_bar))*(h-b2*b2))) {
    //r_bar is too large
    return 1;
  }
  double j = calc_j(b2, r_bar, phi_bar, z);
  double A = calc_A(j, h, b1, b2);
  if (gsl_isnan(A)) {
    cout << "A is nan!" << endl;
  }
  return_values[0] = acos((A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[0], z)) {
    return_values[0] = 10; // 10 means that the theta_bar value is invalid and should be skipped
  }
  return_values[1] = acos((-A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[1], z)) {
    return_values[1] = 10;
  }
  return_values[2] = -acos((A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[2], z)) {
    return_values[2] = 10;
  }
  return_values[3] = -acos((-A+j)/(2*h));
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[3], z)) {
    return_values[3] = 10;
  }
  return 0;
}

double calc_b_bar(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar, double z) {
  return 1/cos(theta_bar)*(b_min + (1-z)*r*cos(phi) - (1-z)*r_bar*cos(theta_bar+phi_bar));
}

double step(double argument) {
  if (argument >= 0) {
    return 1;
  } else {
    return 0;
  }
}

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double W, double Q2, double M_X) {
  double x_pom = (Q2 + M_X*M_X)/(W*W + Q2);
  return get_dipole_amplitude(table, r, b_min, phi, x_pom);
}

double L_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double z, double Q2, double W, double M_X) {
  //static auto t1 = chrono::high_resolution_clock::now();
  if (z*(1-z)*M_X*M_X-m_f*m_f < 0) {
    return 0;
  }
  double total_integrand = 0;
  double theta_bar[4];
  int r_bar_too_large = calc_theta_bar(theta_bar, r, b_min, phi, r_bar, phi_bar, z);
  if (r_bar_too_large == 1) {
    return 0;
  }
  for (int i=0; i<4; i++) {
    if (theta_bar[i] == 10) {
      continue;
    }
    double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i], z);
    double sub_integrand = r*b_min*r_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*M_X*M_X-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))
    *z*(1-z)*4*Q2*z*z*gsl_pow_2(1-z)*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)
    *dipole_amplitude(r, b_min, phi, W, Q2, M_X)*dipole_amplitude(r_bar, b_min_bar, phi_bar, W, Q2, M_X);
    if (gsl_isnan(sub_integrand)) {
      cout << "L sub_integrand " << i << " is nan" << endl;
      /*
      cout << "r=" << setprecision(debug_precision) << r << endl;
      cout << "b_min=" << setprecision(debug_precision) << b_min << endl;
      cout << "phi=" << setprecision(debug_precision) << phi << endl;
      cout << "r_bar=" << setprecision(debug_precision) << r_bar << endl;
      cout << "phi_bar=" << setprecision(debug_precision) << phi_bar << endl;
      cout << "z=" << setprecision(debug_precision) << z << endl;
      cout << "Q2=" << setprecision(debug_precision) << Q2 << endl;
      cout << "x=" << setprecision(debug_precision) << x << endl;
      cout << "beta=" << setprecision(debug_precision) << beta << endl;
      cout << "sub_integrand=" << setprecision(debug_precision) << sub_integrand << endl;
      */
      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  //static auto t2 = chrono::high_resolution_clock::now();
  //auto duration = chrono::duration_cast<chrono::nanoseconds>(t2-t1);
  //cout << duration.count() << endl;

  return total_integrand;
}

double T_integrand(double r, double b_min, double phi, double r_bar, double phi_bar, double z, double Q2, double W, double M_X) {
  //static auto t1 = chrono::high_resolution_clock::now();
  if (z*(1-z)*M_X*M_X-m_f*m_f < 0) {
    return 0;
  }
  double total_integrand = 0;
  double theta_bar[4];
  int r_bar_too_large = calc_theta_bar(theta_bar, r, b_min, phi, r_bar, phi_bar, z);
  if (r_bar_too_large == 1) {
    return 0;
  }
  for (int i=0; i<4; i++) {
    if (theta_bar[i] == 10) {
      continue;
    }
    double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i], z);
    double sub_integrand = r*b_min*r_bar
    *gsl_sf_bessel_J0(sqrt(z*(1-z)*M_X*M_X-m_f*m_f)*sqrt(r*r+r_bar*r_bar-2*r*r_bar*cos(-theta_bar[i]+phi-phi_bar)))*z*(1-z)
    *(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*r_bar)+epsilon2(z, Q2)*(z*z+gsl_pow_2(1-z))*cos(-theta_bar[i]+phi-phi_bar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*r_bar))
    *dipole_amplitude(r, b_min, phi, W, Q2, M_X)*dipole_amplitude(r_bar, b_min_bar, phi_bar, W, Q2, M_X);

    if (gsl_isnan(sub_integrand)) {
      cout << "T sub_integrand " << i << " is nan" << endl;
      /*
      cout << "r=" << setprecision(debug_precision) << r << endl;
      cout << "b_min=" << setprecision(debug_precision) << b_min << endl;
      cout << "phi=" << setprecision(debug_precision) << phi << endl;
      cout << "r_bar=" << setprecision(debug_precision) << r_bar << endl;
      cout << "phi_bar=" << setprecision(debug_precision) << phi_bar << endl;
      cout << "z=" << setprecision(debug_precision) << z << endl;
      cout << "Q2=" << setprecision(debug_precision) << Q2 << endl;
      cout << "x=" << setprecision(debug_precision) << x << endl;
      cout << "beta=" << setprecision(debug_precision) << beta << endl;
      cout << "sub_integrand=" << setprecision(debug_precision) << sub_integrand << endl;
      */
      sub_integrand = 0;
    }

    total_integrand += sub_integrand;
  }
  //static auto t2 = chrono::high_resolution_clock::now();
  //auto duration = chrono::duration_cast<chrono::nanoseconds>(t2-t1);
  //cout << duration.count() << endl;

  return total_integrand;
}

struct parameters {double Q2; double W;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->Q2, par->W, k[6]);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->Q2, par->W, k[6]);
}

struct thread_par_struct
{
  double Q2;
  double W;
  double &sigma;
  double &sigma_error;
  double &sigma_fit;
  thread_par_struct(double a1, double a2, double &a3, double &a4, double &a5) : Q2(a1), W(a2), sigma(a3), sigma_error(a4), sigma_fit(a5) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  const int dim = 7;
  double res, err;

  //double L_integrand(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {
  double xl[dim] = {0, 0, 0, 0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, r_limit, M_PI, 1, 1000};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.W = par.W;
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
    //static auto t1 = chrono::high_resolution_clock::now();
    status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
    if (status != 0) {cout << "integrate_for_L_sigma: " << status << endl; throw (status);}
    //static auto t2 = chrono::high_resolution_clock::now();
    //auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    //cout << "L iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(L_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at W=" << params.W << endl;
  }
  sigma = res;
  sigma_error = err;
  sigma_fit = gsl_monte_vegas_chisq(L_s);
  cout << "L, Q²=" << params.Q2 << ", W=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 7;
  double res, err;

  //double L_integrand(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double z, double Q2, double x, double beta) {
  double xl[dim] = {0, 0, 0, 0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, r_limit, M_PI, 1, 1000};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.W = par.W;
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
    //static auto t1 = chrono::high_resolution_clock::now();
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
    //static auto t2 = chrono::high_resolution_clock::now();
    //auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    //cout << "T iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(T_s) << ", duration: " << duration.count() << endl;
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at W=" << params.W << endl;
  }
  sigma = res;
  sigma_error = err;
  sigma_fit = gsl_monte_vegas_chisq(T_s);
  cout << "T, Q²=" << params.Q2 << ", W=" << params.W << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  const double Q2 = 0;

  if (nucleus_type == "Pb") {
    double r_limit = 657; // 34.64
    double b_min_limit = 328; // 17.32
  } else {
    double r_limit = 34.64;
    double b_min_limit = 17.32;
  }

  const int W_steps = 10;
  const double W_start = 2e1;
  const double W_stop = 2e4;
  const double W_step = 1.0/(W_steps-1)*log10(W_stop/W_start);

  stringstream r_limit_stream;
  r_limit_stream << fixed << setprecision(0) << r_limit;
  string r_limit_string = r_limit_stream.str();

  TString r_limit_filename_string = r_limit_string;

  stringstream b_limit_stream;
  b_limit_stream << fixed << setprecision(0) << b_min_limit;
  string b_limit_string = b_limit_stream.str();

  TString b_limit_filename_string = b_limit_string;

  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
  load_dipole_amplitudes(table, filename);

  TString title, outfile_name;

  TMultiGraph* T_graphs = new TMultiGraph();

  if (print_r_limit) {
    title = "Transverse exclusive "+dipole_amp_type+" "+nucleus_type+" cross section with r limit: "+r_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else if (print_b_min_limit) {
    title = "Transverse exclusive "+dipole_amp_type+" "+nucleus_type+" cross section with b limit: "+b_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else {
    title = "Transverse exclusive "+dipole_amp_type+" "+nucleus_type+" cross section;W (GeV);cross section (mb)";
  }
  T_graphs->SetTitle(title);

  if (print_r_limit) {
    outfile_name = "data/J_LHC_T_exclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".txt";
  } else if (print_b_min_limit) {
    outfile_name = "data/J_LHC_T_exclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".txt";
  } else {
    outfile_name = "data/J_LHC_T_exclusive_"+dipole_amp_type+"_"+nucleus_type+filename_end+".txt";
  }
  ofstream T_output_file(outfile_name);
  T_output_file << "W;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting T integration" << endl;
  
  double T_W_values[W_steps], T_sigma_values[W_steps], T_W_errors[W_steps], T_sigma_errors[W_steps], T_sigma_fits[W_steps];
  thread T_threads[W_steps];

  for (int i=0; i<W_steps; i++) {
    double W = pow(10, log10(W_start) + i*W_step);
    T_W_values[i] = W;
    T_W_errors[i] = 0;
    thread_par_struct par(Q2, W, T_sigma_values[i], T_sigma_errors[i], T_sigma_fits[i]);
    T_threads[i] = thread(integrate_for_T_sigma, par);
  }

  for (int j=0; j<W_steps; j++) {
    T_threads[j].join();
  }

  TGraphErrors* subgraph = new TGraphErrors(W_steps, T_W_values, T_sigma_values, T_W_errors, T_sigma_errors);
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

  TCanvas* T_sigma_canvas = new TCanvas("J_LHC_T_exclusive_sigma_canvas", "", 1100, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  gPad->SetLogy();

  //T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  if (print_r_limit) {
    outfile_name = "figures/J_LHC_T_exclusive_"+dipole_amp_type+"_"+nucleus_type+"_r_"+r_limit_filename_string+filename_end+".pdf";
  } else if (print_b_min_limit) {
    outfile_name = "figures/J_LHC_T_exclusive_"+dipole_amp_type+"_"+nucleus_type+"_b_"+b_limit_filename_string+filename_end+".pdf";
  } else {
    outfile_name = "figures/J_LHC_T_exclusive_"+dipole_amp_type+"_"+nucleus_type+filename_end+".pdf";
  }
  T_sigma_canvas->Print(outfile_name);

  return 0;
}
