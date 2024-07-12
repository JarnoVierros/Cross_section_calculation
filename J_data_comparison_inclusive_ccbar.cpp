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

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "new_dipole_amp_reader.h"

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

static array<array<array<array<double, 4>, 81>, 30>, 30> table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(m_f*m_f + z*(1-z)*Q2);
}

double dipole_amplitude(double r, double b, double x) {
  return 2*M_PI*b*get_dipole_amplitude(table, r, b, x, false);
}

double L_integrand(double r, double z, double b, double Q2, double x) {
  return 2*M_PI*r*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, b, x);
}

double T_integrand(double r, double z, double b, double Q2, double x) {
  return 2*M_PI*r*(m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, b, x);
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], k[2], par->Q2, par->x);
}

bool array_contains_similar(int attay_size, double array[], double element) {
  for (int i=0; i<attay_size; i++) {
    if (0.97 < abs(array[i]/element) && abs(array[i]/element) < 1.03) {
      return true;
    }
  }
  return false;
}

int main() {

  gsl_set_error_handler_off();

  //const double Q2_selection = 2000;

  //const double integration_radius = 100;
  const int warmup_calls = 10000;
  const int integration_calls = 100000;
  const int integration_iterations = 1;

  const string data_filename = "data/HERA_data.dat";

  string dipamp_filename = "data/dipole_amplitude_with_IP_dependence.csv";
  load_dipole_amplitudes(table, dipamp_filename);

  double chisq = 0;
  int ndf = 0;

  double Q2_selections[12] = {2.5, 5, 7, 12, 18, 32, 60, 120, 200, 350, 650, 2000};
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
    cout << "Finished reading file" << endl;

    const int dim = 3;
    double res, err;

    double xl[3] = {0, 0, 0};
    double xu[3] = {34, 1, 20};

    struct parameters params = {1, 1};

    const gsl_rng_type *T;
    gsl_rng *rng;

    gsl_monte_function L_G = {&L_g, dim, &params};
    gsl_monte_function T_G = {&T_g, dim, &params};

    gsl_rng_env_setup ();
    int status = 0;

    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);

    double measured_x[size(Q2_values)], measured_sigma[size(Q2_values)], measured_x_error[size(Q2_values)], measured_sigma_error[size(Q2_values)], model_sigma[size(Q2_values)];

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
      measured_x[j] = temp_measured_x;
      measured_x_error[j] = 0;
      measured_sigma[j] = measured_sigma_values[j];
      measured_sigma_error[j] = relative_measurement_errors[j]/100*measured_sigma_values[j];

      gsl_monte_vegas_state *L_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, warmup_calls, rng, L_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }
      
      double sigma_L = res;

      gsl_monte_vegas_free(L_s);

      gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      for (int i=0; i<integration_iterations; i++) {
        status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
        if (status != 0) {throw "gsl error";}
      }

      double sigma_T = res;

      gsl_monte_vegas_free(T_s);

      double F_L = Q2_values[j]/(4*M_PI*M_PI*alpha_em)*sigma_L;
      double F_T = Q2_values[j]/(4*M_PI*M_PI*alpha_em)*sigma_T;
      double F_2 = F_L + F_T;
      double sigma_r = F_2 - y_values[j]*y_values[j]/(1+gsl_pow_2(1-y_values[j]))*F_L;

      model_sigma[j] = sigma_r;

    }
    cout << "Measurement points: " << size(Q2_values) << endl;
    TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);

    TMultiGraph* comparison_graph = new TMultiGraph();
    stringstream stream;
    stream << fixed << setprecision(0) << Q2_selection;
    TString title = "Reduced cross section fit for Q^{2}="+stream.str()+";x;#sigma_{r} (mb)";
    comparison_graph->SetTitle(title);

    TGraphErrors* measurement_data = new TGraphErrors(size(Q2_values), measured_x, measured_sigma, measured_x_error, measured_sigma_error);
    measurement_data->SetTitle("Reduced cross section fit");
    comparison_graph->Add(measurement_data, "P");

    TGraph* model_fit = new TGraph(size(Q2_values), measured_x, model_sigma);
    comparison_graph->Add(model_fit, "*");

    comparison_graph->GetXaxis()->SetLimits(1e-5, 1e-1);

    comparison_graph->Draw("A");
    
    float location[4];
    if (Q2_selection < 60) {
      location[0] = 0.65;
      location[1] = 0.7;
      location[2] = 0.9;
      location[3] = 0.9;
    } else {
      location[0] = 0.15;
      location[1] = 0.7;
      location[2] = 0.4;
      location[3] = 0.9;
    }
    TLegend* legend = new TLegend(location[0], location[1], location[2], location[3]);
    legend->AddEntry(measurement_data,"Measurement data");
    legend->AddEntry(model_fit,"Model fit", "P");
    legend->SetTextSize(0.04);
    legend->Draw();

    gPad->SetLogx();

    TString outfile_name = "figures/J_data_comparison_Q2_"+stream.str()+".pdf";
    comparison_canvas->Print(outfile_name);

    gsl_rng_free(rng);

    for (long unsigned int j=0; j<size(Q2_values); j++) {
      //if (x_values[j] > 1e-3) {
      //  continue;
      //}
      double delta = (model_sigma[j] - measured_sigma[j])/(measured_sigma_error[j]);
      chisq += delta*delta;
      ndf++;
    }
  }

  cout << "chisq=" << chisq << ", ndf=" << ndf << ", chisq/ndf=" << chisq/ndf << endl;
  
  return 0;
}
