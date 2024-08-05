#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TText.h"
#include "TMinuit.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "direct_dipole_amp_reader.h"

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double normalization = 16*alpha_em*N_c*e_f*e_f/(2*M_PI);

const double r_limit = 34.64; // 34.64
const double b_min_limit = 17.32; // 17.32

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> current_table;
static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> bk_table;
static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> bfkl_table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(m_f*m_f + z*(1-z)*Q2);
}

double dipole_amplitude(double r, double b, double phi, double x) {
  return get_dipole_amplitude(current_table, r, b, phi, x);
}

double L_integrand(double r, double b, double phi, double z, double Q2, double x) {
  return r*b*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, b, phi, x);
}

double T_integrand(double r, double b, double phi, double z, double Q2, double x) {
  return r*b*(m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, b, phi, x);
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], k[3], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], k[3], par->Q2, par->x);
}
/*
bool array_contains_similar(int attay_size, double array[], double element) {
  for (int i=0; i<attay_size; i++) {
    if (0.97 < abs(array[i]/element) && abs(array[i]/element) < 1.03) {
      return true;
    }
  }
  return false;
}
*/
void fit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  double normalization = par[0];

  gsl_set_error_handler_off();

  //const double Q2_selection = 2000;

  //const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const string data_filename = "data/HERA_data.dat";

  string bk_dipamp_filename = "data/dipole_amplitude_with_IP_dependence.csv";
  load_dipole_amplitudes(bk_table, bk_dipamp_filename);

  string bfkl_dipamp_filename = "data/dipole_amplitude_with_IP_dependence_bfkl.csv";
  load_dipole_amplitudes(bfkl_table, bfkl_dipamp_filename);

  current_table = bk_table;

  double chisq = 0;
  int ndf = 0;

  double Q2_selections[12] = {2.5, 5, 7, 12, 18, 32, 60, 120, 200, 350, 650, 2000};
  
  TMultiGraph* comparison_graphs[size(Q2_selections)];
  TGraphErrors* measurement_datas[size(Q2_selections)];
  TGraph* bk_model_fits[size(Q2_selections)];
  TGraph* bfkl_model_fits[size(Q2_selections)];
  
  for (int n=0; n<12; n++) {

    double Q2_selection = Q2_selections[n];

    vector<double> Q2_values;
    vector<double> x_values;
    vector<double> y_values;
    vector<double> measured_sigma_values;
    vector<double> relative_measurement_errors;

    ifstream data_file(data_filename);

    cout << "Reading: " << data_filename << endl;
    string line;
    while(getline (data_file, line)) {
      long unsigned int i = 0;
      string value = "";
      while(line[i] != ' ') {
        value += line[i];
        i++;
      }
      if (stod(value) != Q2_selection) {continue;}
      Q2_values.push_back(stod(value));
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
      y_values.push_back(stod(value));
      i++;

      value = "";
      while(line[i] != ' ') {
        value += line[i];
        i++;
      }
      measured_sigma_values.push_back(stod(value));
      i++;

      for (int j=0; j<3; j++) {
        while(line[i] != ' ') {
          i++;
        }
        i++;
      }

      value = "";
      while(i<line.length()) {
        value += line[i];
        i++;
      }
      relative_measurement_errors.push_back(stod(value));
    }

    if (n == 11) {
      Q2_values.push_back(2000);
      x_values.push_back(0.04000);
      y_values.push_back(0.395);
      measured_sigma_values.push_back(-1);
      relative_measurement_errors.push_back(0);

      Q2_values.push_back(2000);
      x_values.push_back(0.07000);
      y_values.push_back(0.395);
      measured_sigma_values.push_back(-1);
      relative_measurement_errors.push_back(0);
    }

    cout << "Finished reading file" << endl;

    const int dim = 4;
    double res, err;

    double xl[4] = {0, 0, 0, 0};
    double xu[4] = {r_limit, b_min_limit, M_PI, 1};

    struct parameters params = {1, 1};

    const gsl_rng_type *T;
    gsl_rng *rng;

    gsl_monte_function L_G = {&L_g, dim, &params};
    gsl_monte_function T_G = {&T_g, dim, &params};

    gsl_rng_env_setup ();
    int status = 0;

    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);

    double measured_x[size(Q2_values)], measured_sigma[size(Q2_values)], measured_x_error[size(Q2_values)], measured_sigma_error[size(Q2_values)], bk_model_sigma[size(Q2_values)], bfkl_model_sigma[size(Q2_values)];

    for (long unsigned int j=0; j<size(Q2_values); j++) {
      /*
      if (x_values[j] != 0.00080) {
        measured_x[j] = x_values[j];
        measured_x_error[j] = 0;
        measured_sigma[j] = 0;
        measured_sigma_error[j] = 0;
        model_sigma[j] = 0;
        continue;
      }
      */
      cout << "Integrating at Q^2=" << Q2_values[j] << endl;

      params.Q2 = Q2_values[j];
      params.x = x_values[j];
      /*
      int depth = 0;

      double temp_measured_x = x_values[j];
      double multiplier = 1;
      
      while (array_contains_similar(size(Q2_values), measured_x, temp_measured_x)) {
        cout << "x=" << temp_measured_x << " already occupied" << endl;
        multiplier += 0.1;
        temp_measured_x = x_values[j]*multiplier;
        cout << "changed to x=" << temp_measured_x << endl;
        depth++;
        if (depth>100) {return 1;}
      }
      */
      measured_x[j] = x_values[j];
      measured_x_error[j] = 0;
      measured_sigma[j] = measured_sigma_values[j];
      measured_sigma_error[j] = relative_measurement_errors[j]/100*measured_sigma_values[j];

      ///BK model

      current_table = bk_table;

      gsl_monte_vegas_state *bk_L_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, bk_L_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, bk_L_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }
      
      double sigma_L = normalization*res;

      gsl_monte_vegas_free(bk_L_s);

      gsl_monte_vegas_state *bk_T_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, bk_T_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, bk_T_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }

      double sigma_T = normalization*res;

      gsl_monte_vegas_free(bk_T_s);

      double F_L = Q2_values[j]/(4*M_PI*M_PI*alpha_em)*sigma_L;
      double F_T = Q2_values[j]/(4*M_PI*M_PI*alpha_em)*sigma_T;
      double F_2 = F_L + F_T;
      double sigma_r = F_2 - y_values[j]*y_values[j]/(1+gsl_pow_2(1-y_values[j]))*F_L;

      bk_model_sigma[j] = sigma_r;

    }

    gsl_rng_free(rng);

    for (long unsigned int j=0; j<size(Q2_values); j++) {
      if (measured_x[j] < 1e-2 && Q2_selection <= 60) {
        double delta = (bk_model_sigma[j] - measured_sigma[j])/(measured_sigma_error[j]);
        chisq += delta*delta;
        ndf++;
      }
    }
  }

  cout << "chisq=" << chisq << ", ndf=" << ndf << ", chisq/ndf=" << chisq/ndf << endl;
  
  f = chisq;
}


int main() {

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  Double_t val1,err1;
  Double_t arglist[10];
  Int_t ierflg = 0;
  static Double_t vstart[1];
  static Double_t step[1];

  TMinuit* gMinuit = new TMinuit(3);
  
  gMinuit->SetFCN(fit);

  arglist[0] = 1;

  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  vstart[0] = 1;

  step[0] =0.1;
  gMinuit->mnparm(0, "a0", vstart[0], step[0], 0,0,ierflg);

  arglist[0] = 500;
  arglist[1] = 1.;
  cout << "Starting fitting" << endl;
  //gMinuit->mnhess();
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  //gMinuit->mnhess();

  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  gMinuit->GetParameter(0, val1, err1);

  gMinuit ->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  cout << "val1: " << val1 << ", err1: " << err1 << endl;
  cout << "chi2: " << amin << endl;
  
}