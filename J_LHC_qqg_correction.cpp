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
#include "Jani_qqg_contribution_calculator.cpp"

const double alpha_em = 1.0/137;
const int N_c = 3;
//const double e_f = 1.0/3; //2.0/3
//const double m_f = 4.18; //1.27 GeV
const double e_f = 2.0/3; //
const double m_f = 1.27; // GeV

const int warmup_calls = 10000;
const int integration_calls = 100000;
const int integration_iterations = 1;

const string dipole_amp_type = "bk";
const string nucleus_type = "Pb";
const string diffraction = "_diffraction";//_diffraction this tells it to use diffractive dipole amplitude
const string filename_end = "";
const string particle_name = "c";

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

static InterpMultilinear<4, double>* interpolator;

struct parameters {double Q2; double W;};

struct thread_par_struct
{
  double Q2;
  double W;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double a2, double &a3, double &a4) : Q2(a1), W(a2), sigma(a3), sigma_error(a4) {}
};

struct integrand_parameters {
  double W;
};

double integrand(double beta, void * parameters) {
  struct integrand_parameters * params = (struct integrand_parameters *)parameters;
  double W = params->W;

  double xpom = 4*m_f*m_f/(beta*W*W);

  double result, error, fit;

  calc_total_xpomF_Tqqg_contribution(beta, xpom, 0, result, error, fit);

  return result;
}

double integrate_qqg_contribution(thread_par_struct par) {

  double W = par.W;

  gsl_integration_workspace * wspace = gsl_integration_workspace_alloc (10000);

  struct integrand_parameters parameters = {W};

  double result, error;

  gsl_function F;
  F.function = &integrand;
  F.params = &parameters;

  gsl_integration_qag(&F, 0, 1, 0, 0.01, 1000, 4, wspace, &result, &error);

  par.sigma = result;
  par.sigma_error = error;

  gsl_integration_workspace_free(wspace);

  //cout << "Result: " << result << ", error: " << error << endl;
  return result;
}

int main() {

  gsl_set_error_handler_off();

  if (nucleus_type == "Pb") {
    r_limit = 657; // 34.64
    b_min_limit = 328; // 17.32
  } else if (nucleus_type == "p") {
    r_limit = 34.64;
    b_min_limit = 17.32;
  } else {
    cout << "invalid nucleus type" << endl;
    throw 1;
  }

  const int Q2_values[] = {0};

  
  const int W_steps = 50;
  const double W_start = 3e1;
  const double W_stop = 2e4;
  const double W_step = 1.0/(W_steps-1)*log10(W_stop/W_start);

  string filename;
  if (dipole_amp_type == "vector") {
    filename = "data/bk_p_mu02_0.66.csv";
  } else {
    filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+diffraction+".csv";
  }
  if (nucleus_type == "p") {
    if (diffraction=="_diffraction" && dipole_amp_type == "bfkl") {
      load_Pb_dipole_amplitudes(Pb_table, filename);
      create_Pb_interpolator(Pb_table, interpolator);
    } else {
      load_p_dipole_amplitudes(p_table, filename);
      create_p_interpolator(p_table, interpolator);
    }
  } else if (nucleus_type == "Pb") {
    load_Pb_dipole_amplitudes(Pb_table, filename);
    create_Pb_interpolator(Pb_table, interpolator);
  } else {
    throw 1;
  }

  TMultiGraph* T_graphs = new TMultiGraph();
  T_graphs->SetTitle("Diffractive transverse cross section;W (GeV);cross section (mb)");

  ofstream T_output_file("data/LHC_qqg_correction_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+".txt");
  T_output_file << "W (GeV);sigma (mb);sigma error (mb)" << endl;

  cout << "Starting T integration" << endl;
  
  for (long unsigned int j=0; j<size(Q2_values); j++) {
    double T_W_values[W_steps], T_sigma_values[W_steps], T_W_errors[W_steps], T_sigma_errors[W_steps];
    thread T_threads[W_steps];

    for (int i=0; i<W_steps; i++) {
      double W = pow(10, log10(W_start) + i*W_step);
      T_W_values[i] = W;
      T_W_errors[i] = 0;
      thread_par_struct par(Q2_values[j], W, T_sigma_values[i], T_sigma_errors[i]);
      T_threads[i] = thread(integrate_qqg_contribution, par);
      //this_thread::sleep_for(30s);
    }

    for (int j=0; j<W_steps; j++) {
      T_threads[j].join();
    }

    TGraphErrors* subgraph = new TGraphErrors(W_steps, T_W_values, T_sigma_values, T_W_errors, T_sigma_errors);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[j]);
    subgraph->SetTitle(subgraph_name);
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
  }
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);
  T_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  gPad->SetLogy();

  T_sigma_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);

  TString fig_filename = "figures/LHC_qqg_correction_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+".pdf";
  T_sigma_canvas->Print(fig_filename);

  return 0;
}
