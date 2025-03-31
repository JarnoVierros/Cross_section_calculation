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
#include "TLatex.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <thread>
using namespace std;

#include "linterp.h"
#include "direct_dipole_amp_reader.h"

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3; //2.0/3
const double m_f = 1.27; //1.27 GeV

const double m_Pb = (207.2-82*5.486e-4)*931.49410372*1.0/1000;
const double p_Pb = sqrt(pow(5.36*1000, 2)/4 - m_Pb*m_Pb);

static const double global_Q2 = 0;

const double normalization = 2/pow(2*M_PI, 4)*alpha_em*N_c*e_f*e_f;

static double r_limit; // 34.64
static double b_min_limit; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;
const string dipole_amp_type = "bk";
const string nucleus_type = "Pb";
const string diffraction = "";//_diffraction
const string filename_end = "_converging";//_1mil
const string particle_name = "c";

const int warmup_calls = 100000;
const int integration_calls = 100000000;
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

double shifted_x_at_Q2_0(double z, double pT, double y, double Q2) {
  if (Q2 != 0) {
    cout << "the shifted x formula being used applies only ato Q2=0" << endl;
    throw 1;
  }

  double p_gamma = m_f/(z*sqrt(2))*(sqrt(exp(2*y)+1) + sqrt(exp(2*y)-1-2*pT*pT/(m_f*m_f)));
  double shifted_x = 2*m_f*m_f/(p_gamma*(sqrt(m_Pb*m_Pb+p_Pb*p_Pb) - p_Pb));
  return shifted_x;
}

double dipole_amplitude(double z, double r, double b_min, double phi, double pT, double y, double Q2) {
  double shifted_x = shifted_x_at_Q2_0(z, pT, y, Q2);
  //cout << "x: " << shifted_x << endl;

  if (calc_max_phi(r, b_min) < phi) {
    return 0;
  } else {
    array<double, 4> args = {log(r), log(b_min), phi, log(shifted_x)};
    return exp(interpolator->interp(args.begin()));
  }
}

double dipole_amplitude_from_xy(double rx, double Tx, double ry, double Ty, double z, double pT, double rapidity, double Q2) {
  double dip_r, dip_bmin, dip_phi;
  if (rx < ry) {
    dip_r = sqrt(rx*rx-2*rx*ry*cos(Tx-Ty)+ry*ry);
    dip_bmin = rx;
    dip_phi = acos((rx*cos(Tx)*(ry*cos(Ty)-rx*cos(Tx))+rx*sin(Tx)*(ry*sin(Ty)-rx*sin(Tx)))/(rx*dip_r));
  } else {
    dip_r = sqrt(ry*ry-2*ry*rx*cos(Ty-Tx)+rx*rx);
    dip_bmin = ry;
    dip_phi = acos((ry*cos(Ty)*(rx*cos(Tx)-ry*cos(Ty))+ry*sin(Ty)*(rx*sin(Tx)-ry*sin(Ty)))/(ry*dip_r));
  }
  return dipole_amplitude(z, dip_r, dip_bmin, dip_phi, pT, rapidity, Q2);
}

double T_integrand(double rx, double Tx, double ry, double Ty, double rxbar, double Txbar, double z, double Q2, double pT, double rapidity) {
  double r = sqrt(rx*rx-2*rx*ry*cos(Tx-Ty)+ry*ry);
  double rbar = sqrt(rxbar*rxbar-2*rxbar*ry*cos(Txbar-Ty)+ry*ry);
  double dip_amp_1, dip_amp_2, dip_amp_3;
  dip_amp_1 = dipole_amplitude_from_xy(rx, Tx, ry, Ty, z, pT, rapidity, Q2);
  dip_amp_2 = dipole_amplitude_from_xy(ry, Ty, rxbar, Txbar, z, pT, rapidity, Q2);
  dip_amp_3 = dipole_amplitude_from_xy(rx, Tx, rxbar, Txbar, z, pT, rapidity, Q2);
  return rx*ry*rxbar*cos(pT*(rx*cos(Tx)-rxbar*cos(Txbar)))*(m_f*m_f*gsl_sf_bessel_K0(epsilon(z, Q2)*r)*gsl_sf_bessel_K0(epsilon(z, Q2)*rbar) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*(rx*rxbar*cos(Tx-Txbar)-rx*ry*cos(Tx-Ty)-ry*rxbar*cos(Ty-Txbar)+ry*ry)/(r*rbar)*gsl_sf_bessel_K1(epsilon(z, Q2)*r)*gsl_sf_bessel_K1(epsilon(z, Q2)*rbar))*(dip_amp_1+dip_amp_2-dip_amp_3);
}

struct parameters {double pT; double y;};

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], k[4], k[5], k[6], global_Q2, par->pT, par->y);
}

struct thread_par_struct
{
  double pT;
  double y;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double a2, double &a3, double &a4) : pT(a1), y(a2), sigma(a3), sigma_error(a4) {}
};

void integrate_for_T_sigma(thread_par_struct par) {

  const int dim = 7;
  double res, err;

  const double integration_radius = 200;
  double xl[dim] = {0, 0, 0, 0, 0, 0, 0};
  double xu[dim] = {integration_radius, 2*M_PI, integration_radius, 2*M_PI, integration_radius, 2*M_PI, 1};

  struct parameters params = {1, 1};
  params.pT = par.pT;
  params.y = par.y;
  double &sigma = par.sigma;
  double &sigma_error = par.sigma_error;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function T_G = {&T_g, dim, &params};

  gsl_rng_env_setup ();
  int status = 0;

  //T = gsl_rng_default;
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
  gsl_rng_set(rng, 246187246812);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  for (int i=0; i<integration_iterations; i++) {
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
    if (status != 0) {cout << "integrate_for_T_sigma: " << status << endl; throw (status);}
  }
  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at pT=" << params.pT << ", y=" << params.y << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "T, pT=" << params.pT << ", y=" << params.y << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+diffraction+".csv";
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

  const int y_steps = 7;
  const double y_start = 0;
  const double y_stop = 20;
  const double y_step = (y_stop-y_start)/y_steps;

  const int pT_steps = 1;
  double pT_values[pT_steps] = {0};
  /*
  const double pT_start = 1e0;
  const double pT_stop = 1e3;
  const double pT_step = 1.0/(pT_steps-1)*log10(pT_stop/pT_start);
  */

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

  TMultiGraph* T_graphs = new TMultiGraph();

  TString title, outfile_name;
  if (print_r_limit) {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+diffraction+" cross section with r limit: "+r_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else if (print_b_min_limit) {
    title = "Transverse inclusive "+dipole_amp_type+" "+nucleus_type+diffraction+" cross section with b limit: "+b_limit_string+" GeV^-1;W (GeV);cross section (mb)";
  } else {
    title = "Inclusive D0 production in #gamma "+nucleus_type+diffraction+" scattering using "+dipole_amp_type+" evolution;y;cross section (GeV^-2)";
  }
  T_graphs->SetTitle(title);

  if (print_r_limit) {
    outfile_name = "output/J_LHC_T_inclusive_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+"_r_"+r_limit_filename_string+filename_end+".txt";
  } else if (print_b_min_limit) {
    outfile_name = "output/J_LHC_T_inclusive_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+"_b_"+b_limit_filename_string+filename_end+".txt";
  } else {
    outfile_name = "output/rapidity_LHC_inclusive_D0_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+filename_end+".txt";
  }
  ofstream T_output_file(outfile_name);
  T_output_file << "pT;y;cross section (GeV^-2);cross section error (GeV^-2)" << endl;

  cout << "Starting T integration" << endl;
  
  double T_pT_values[pT_steps*y_steps], T_y_values[pT_steps*y_steps], T_sigma_values[pT_steps*y_steps], T_pT_errors[pT_steps*y_steps], T_y_errors[pT_steps*y_steps], T_sigma_errors[pT_steps*y_steps];
  thread T_threads[pT_steps*y_steps];

  const int max_threads = 7;
  int active_threads = 0;
  for (int i=0; i<pT_steps; i++) {
    for (int j=0; j<y_steps; j++) {
      //double pT = pow(10, log10(pT_start) + i*pT_step);
      //cout << "pT=" << pT << endl;
      double pT = pT_values[i];
      double y = y_start + j*y_step;
      T_pT_values[i*y_steps+j] = pT;
      T_pT_errors[i*y_steps+j] = 0;
      T_y_values[i*y_steps+j] = y;
      T_y_errors[i*y_steps+j] = 0;
      thread_par_struct par(pT, y, T_sigma_values[i*y_steps+j], T_sigma_errors[i*y_steps+j]);
      T_threads[i*y_steps+j] = thread(integrate_for_T_sigma, par);
      active_threads += 1;

      if (active_threads == max_threads) {
        for (int k=0; k<active_threads; k++) {
          T_threads[i*y_steps+j-k].join();
        }
        active_threads = 0;
      }
    }
  }
  for (int k=0; k<active_threads; k++) {
    T_threads[pT_steps*y_steps-k-1].join();
  }
  /*
  for (int j=0; j<pT_steps*y_steps; j++) {
    T_threads[j].join();
  }
  */

  for (int i=0; i<pT_steps; i++) {

    double y[y_steps], sigma[y_steps], y_error[y_steps], sigma_error[y_steps];
    for (int j=0; j<y_steps; j++) {
      y[j] = T_y_values[i*y_steps+j];
      sigma[j] = T_sigma_values[i*y_steps+j];
      y_error[j] = T_y_errors[i*y_steps+j];
      sigma_error[j] = T_sigma_errors[i*y_steps+j];
    }

    TGraphErrors* subgraph = new TGraphErrors(y_steps, y, sigma, y_error, sigma_error);
    TString subgraph_title = "pT="+to_string(pT_values[i]);
    subgraph->SetTitle(subgraph_title);
    T_graphs->Add(subgraph);
  }

  for (int i=0; i<pT_steps*y_steps; i++) {
    ostringstream pT;
    pT << T_pT_values[i];
    ostringstream y;
    y << T_y_values[i];
    ostringstream sigma;
    sigma << T_sigma_values[i];
    ostringstream sigma_err;
    sigma_err << T_sigma_errors[i];
    string line = pT.str() + ";" + y.str() + ";" + sigma.str() + ";" + sigma_err.str();
    T_output_file << line << endl;
  }
  
  T_output_file.close();

  TCanvas* T_sigma_canvas = new TCanvas("rapidity_LHC_inclusive_D0_canvas", "", 1100, 600);
  T_graphs->Draw("A PMC PLC");

  //gPad->SetLogx();
  //gPad->SetLogy();

  if (size(pT_values) > 1) {
    T_sigma_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  }

  //TLatex* Q2_text = new TLatex(1.5e0, 2e3, "Q^{2} = 0 GeV^{2}");
  //Q2_text->Draw("same");

  if (print_r_limit) {
    outfile_name = "figures/rapidity_LHC_inclusive_D0_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+"_r_"+r_limit_filename_string+filename_end+".pdf";
  } else if (print_b_min_limit) {
    outfile_name = "figures/rapidity_LHC_inclusive_D0_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+"_b_"+b_limit_filename_string+filename_end+".pdf";
  } else {
    outfile_name = "figures/rapidity_LHC_inclusive_D0_"+particle_name+"_"+dipole_amp_type+"_"+nucleus_type+diffraction+filename_end+".pdf";
  }
  T_sigma_canvas->Print(outfile_name);

  return 0;
}