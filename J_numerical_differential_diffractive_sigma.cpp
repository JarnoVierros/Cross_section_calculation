
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
const double m_f = 1.27; //GeV 1.27

const double normalization = 16/gsl_pow_3(2*M_PI)*alpha_em*N_c*e_f*e_f;

const double r_limit = 34.64; // 34.64
const double b_min_limit = 17.32; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;

const int warmup_calls = 10000;
const int integration_calls = 100000;
const int integration_iterations = 1;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

void handle_status(int status) {
  if (status != 0) {
    if (status == 18) {
      cout << "status 18" << endl;
    } else if (status == 22) {
      cout << "warning: divergent qags integral" << endl;
    } else {
      cout << "L_integrand: " << status << endl;
      throw (status);
    }
  }
}

double calc_h(double r, double b_min, double phi) {
  return 4*b_min*b_min + 4*b_min*r*cos(phi) + r*r;
}

double calc_b1(double r, double b_min, double phi, double theta) {
  return b_min*cos(theta) + r/2*cos(theta+phi);
}

double calc_b2(double r, double b_min, double phi, double theta) {
  return b_min*sin(theta) + r/2*sin(theta+phi);
}

double calc_j(double b2, double r_bar, double phi_bar) {
  return 4*b2*r_bar*sin(phi_bar);
}

double calc_A(double j, double h, double b1, double b2) {
  return sqrt(j*j + h*(16*b1*b1-j*j/(b2*b2)));
}

double calc_theta_bar(double return_values[4], double r, double b_min, double phi, double theta, double r_bar, double phi_bar) {
  double b1 = calc_b1(r, b_min, phi, theta);
  double b2 = calc_b2(r, b_min, phi, theta);
  double h = calc_h(r, b_min, phi);
  double j = calc_j(b2, r_bar, phi_bar);
  double A = calc_A(j, h, b1, b2);
  return_values[0] = acos((A+j)/(2*h));
  return_values[1] = acos((-A+j)/(2*h));
  return_values[2] = -acos((A+j)/(2*h));
  return_values[3] = -acos((-A+j)/(2*h));
}

double calc_b_bar(double r, double b_min, double phi, double theta, double r_bar, double phi_bar, double theta_bar) {
  return 1/cos(theta_bar)*(b_min*cos(theta) + r/2*cos(theta+phi) - r_bar/2*cos(theta_bar+phi_bar));
}


double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double b_min, double phi, double x) {
  return get_dipole_amplitude(table, r, b_min, phi, x);
}

struct integration_parameters {double r; double b_min; double phi; double theta; double r_bar; double phi_bar; double z; double Q2; double x_pom; double beta;};

double phibar_integrand(double phibar, void* phibar_params) {
  struct integration_parameters *par = (struct integration_parameters *)phibar_params;

  double total_integrand = 0;
  double theta_bar[4];
  calc_theta_bar(theta_bar, par->r, par->b_min, par->phi, par->theta, par->r_bar, par->phi_bar);
  for (int i=0; i<4; i++) {
    double b_min_bar = calc_b_bar(par->r, par->b_min, par->phi, par->theta, par->r_bar, par->phi_bar, theta_bar[i]);
    total_integrand += dipole_amplitude(par->r_bar, b_min_bar, par->phi_bar, par->beta*par->x_pom)
    *gsl_sf_bessel_J0(sqrt(par->z*(1-par->z)*par->Q2*(1/par->beta-1)-m_f*m_f)*sqrt(gsl_pow_2(par->r*cos(par->phi+par->theta)-par->r_bar*cos(theta_bar[i]+par->phi_bar)) + gsl_pow_2(par->r*sin(par->phi+par->theta)-par->r_bar*sin(par->phi_bar+theta_bar[i]))))
    *4*par->Q2*par->z*par->z*gsl_pow_2(1-par->z)*gsl_sf_bessel_K0(epsilon(par->z, par->Q2)*par->r)*gsl_sf_bessel_K0(epsilon(par->z, par->Q2)*par->r_bar);
  }
  return total_integrand;
}

double rbar_integrand(double rbar, void* rbar_params) {
  struct integration_parameters *par = (struct integration_parameters *)rbar_params;
  gsl_integration_workspace * w_phibar = gsl_integration_workspace_alloc (1000);
  double result, error;
  struct integration_parameters phibar_params = {par->r, par->b_min, par->phi, par->theta, rbar, 0, par->z, par->Q2, par->x_pom, par->beta};
  gsl_function F_phibar;
  F_phibar.function = &phibar_integrand;
  F_phibar.params = &phibar_params;

  int status = gsl_integration_qags (&F_phibar, 0, r_limit, 0, 0.01, 1000, w_phibar, &result, &error);
  handle_status(status);

  gsl_integration_workspace_free(w_phibar);
  return rbar*result;
}

double theta_integrand(double theta, void* theta_params) {
  struct integration_parameters *par = (struct integration_parameters *)theta_params;
  gsl_integration_workspace * w_rbar = gsl_integration_workspace_alloc (1000);
  double result, error;
  struct integration_parameters rbar_params = {par->r, par->b_min, par->phi, theta, 0, 0, par->z, par->Q2, par->x_pom, par->beta};
  gsl_function F_rbar;
  F_rbar.function = &rbar_integrand;
  F_rbar.params = &rbar_params;

  int status = gsl_integration_qags (&F_rbar, 0, 2*M_PI, 0, 0.01, 1000, w_rbar, &result, &error);
  handle_status(status);

  gsl_integration_workspace_free(w_rbar);
  return dipole_amplitude(par->r, par->b_min, par->phi, par->beta*par->x_pom)*result;
}

double phi_integrand(double phi, void* phi_params) {
  struct integration_parameters *par = (struct integration_parameters *)phi_params;
  gsl_integration_workspace * w_theta = gsl_integration_workspace_alloc (1000);
  double result, error;
  struct integration_parameters theta_params = {par->r, par->b_min, phi, 0, 0, 0, par->z, par->Q2, par->x_pom, par->beta};
  gsl_function F_theta;
  F_theta.function = &theta_integrand;
  F_theta.params = &theta_params;

  int status = gsl_integration_qags (&F_theta, 0, M_PI, 0, 0.01, 1000, w_theta, &result, &error);
  handle_status(status);

  gsl_integration_workspace_free(w_theta);
  return result;
}

double bmin_integrand(double bmin, void* bmin_params) {
  struct integration_parameters *par = (struct integration_parameters *)bmin_params;
  gsl_integration_workspace * w_phi = gsl_integration_workspace_alloc (1000);
  double result, error;
  struct integration_parameters phi_params = {par->r, bmin, 0, 0, 0, 0, par->z, par->Q2, par->x_pom, par->beta};
  gsl_function F_phi;
  F_phi.function = &phi_integrand;
  F_phi.params = &phi_params;

  int status = gsl_integration_qags (&F_phi, 0, M_PI, 0, 0.01, 1000, w_phi, &result, &error);
  handle_status(status);

  gsl_integration_workspace_free(w_phi);
  return bmin*result;
}

double r_integrand(double r, void* r_params) {
  struct integration_parameters *par = (struct integration_parameters *)r_params;
  gsl_integration_workspace * w_bmin = gsl_integration_workspace_alloc (1000);
  double result, error;
  struct integration_parameters bmin_params = {r, 0, 0, 0, 0, 0, par->z, par->Q2, par->x_pom, par->beta};
  gsl_function F_bmin;
  F_bmin.function = &bmin_integrand;
  F_bmin.params = &bmin_params;

  int status = gsl_integration_qags (&F_bmin, 0, b_min_limit, 0, 0.01, 1000, w_bmin, &result, &error);
  handle_status(status);

  gsl_integration_workspace_free(w_bmin);
  return r*result;
}

double z_integrand(double z, void* z_params) {
  struct integration_parameters *par = (struct integration_parameters *)z_params;
  
  if (z*(1-z)*par->Q2*(1/par->beta-1)-m_f*m_f < 0) {
    return 0;
  }
  
  gsl_integration_workspace * w_r = gsl_integration_workspace_alloc (1000);
  double result, error;
  struct integration_parameters r_params = {0, 0, 0, 0, 0, 0, z, par->Q2, par->x_pom, par->beta};
  gsl_function F_r;
  F_r.function = &r_integrand;
  F_r.params = &r_params;

  int status = gsl_integration_qags (&F_r, 0, r_limit, 0, 0.01, 1000, w_r, &result, &error);
  handle_status(status);

  gsl_integration_workspace_free(w_r);
  return z*(1-z)*result;
}

void L_integral(double &sigma, double &sigma_error, double Q2, double x_pom, double beta) {
  gsl_integration_workspace * w_z = gsl_integration_workspace_alloc (1000);
  double result, error;
  struct integration_parameters z_params = {0, 0, 0, 0, 0, 0, 0, Q2, x_pom, beta};
  gsl_function F_z;
  F_z.function = &z_integrand;
  F_z.params = &z_params;

  int status = gsl_integration_qags (&F_z, 0, 1, 0, 0.01, 1000, w_z, &result, &error);
  handle_status(status);

  gsl_integration_workspace_free(w_z);
  
  sigma = result;
  sigma_error = error;
}

int main() {

  gsl_set_error_handler_off();

  double sigma;
  double sigma_error;
  L_integral(sigma, sigma_error, 2.5, 0.01, 0.1);
  cout << "sigma: " << sigma << ", error: " << sigma_error << endl;

  /*

  const int Q2_values[] = {1, 3, 5, 8, 10};

  const int x_steps = 30;
  const double x_start = 1e-5;
  const double x_stop = 0.01;
  const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

  string filename = "data/dipole_amplitude_with_IP_dependence.csv";
  load_dipole_amplitudes(table, filename);

  TMultiGraph* L_graphs = new TMultiGraph();
  L_graphs->SetTitle("Diffractive longitudinal cross section;x;cross section (mb)");

  ofstream L_output_file("data/diff_L_sigma_x.txt");
  L_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting L integration" << endl;
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double L_x_values[x_steps], L_sigma_values[x_steps], L_x_errors[x_steps], L_sigma_errors[x_steps];
    thread L_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      L_x_values[i] = x;
      L_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, L_sigma_values[i], L_sigma_errors[i]);
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
      string line = to_string(Q2_values[j]) + ";" + x.str() + ";" + sigma.str() + ";" + sigma_err.str();
      L_output_file << line << endl;
    }
  }
  L_output_file.close();

  TCanvas* L_sigma_canvas = new TCanvas("diff_L_sigma_canvas", "", 1000, 600);
  L_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  L_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  L_sigma_canvas->Print("figures/diff_L_sigma_x_distribution.pdf");
 

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Diffractive transverse cross section;x;cross section (mb)");

  ofstream T_output_file("data/diff_T_sigma_x.txt");
  T_output_file << "Q2 (GeV);x;sigma (mb);sigma error (mb)" << endl;

  cout << "Starting T integration" << endl;
  
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_x_values[x_steps], T_sigma_values[x_steps], T_x_errors[x_steps], T_sigma_errors[x_steps];
    thread T_threads[x_steps];

    for (int i=0; i<x_steps; i++) {
      double x = pow(10, log10(x_start) + i*x_step);
      T_x_values[i] = x;
      T_x_errors[i] = 0;
      thread_par_struct par(Q2_values[j], x, T_sigma_values[i], T_sigma_errors[i]);
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

  T_sigma_canvas->Print("figures/diff_T_sigma_x_distribution.pdf");

  return 0;
  */
}
