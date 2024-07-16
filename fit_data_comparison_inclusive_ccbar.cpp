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


const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV

const double sigma_0 = 2.61843e+01; //2.99416e+01; //mb
const double Q_0 = 1; //GeV
const double x_0 = 8.66918e-05; //7.67079e-05;
const double lambda_star = 3.15655e-01; //3.64361e-01;

const double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double dipole_amplitude(double r, double x) {
  return sigma_0*(1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow(x/x_0, lambda_star/2)))));
}

double L_integrand(double r, double z, double Q2, double x) {
  int status = 0;
  gsl_sf_result bessel_result;
  status = gsl_sf_bessel_K0_e(epsilon(z, Q2)*r, &bessel_result);
  if (status !=0) {
    if (status == 15) {
      bessel_result.val = 0;
    } else {
      cout << "GSL error in L_integrand: " << status << endl;
      throw 1;
    }
  }
  return 2*M_PI*r*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(bessel_result.val)*dipole_amplitude(r, x);
}

double T_integrand(double r, double z, double Q2, double x) {
  int status = 0;
  gsl_sf_result bessel_result_1;
  gsl_sf_result bessel_result_2;
  status = gsl_sf_bessel_K0_e(epsilon(z, Q2)*r, &bessel_result_1);
  if (status !=0) {
    if (status == 15) {
      bessel_result_1.val = 0;
    } else {
      cout << "GSL error in L_integrand: " << status << endl;
      throw 1;
    }
  }
  status = gsl_sf_bessel_K1_e(epsilon(z, Q2)*r, &bessel_result_2);
  if (status !=0) {
    if (status == 15) {
      bessel_result_2.val = 0;
    } else {
      cout << "GSL error in L_integrand: " << status << endl;
      throw 1;
    }
  }
  return 2*M_PI*r*(m_f*m_f*gsl_pow_2(bessel_result_1.val) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*gsl_pow_2(bessel_result_2.val))*dipole_amplitude(r, x);
}

struct parameters {double Q2; double x;};

double L_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*L_integrand(k[0], k[1], par->Q2, par->x);
}

double T_g(double *k, size_t dim, void * params) {
  struct parameters *par = (struct parameters *)params;
  return normalization*T_integrand(k[0], k[1], par->Q2, par->x);
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

  const int warmup_calls = 100000;
  const int integration_calls = 1000000;
  const int integration_iterations = 1;

  const string filename = "data/HERA_data.dat";

  vector<double> Q2_values;
  vector<double> x_values;
  vector<double> y_values;
  vector<double> measured_sigma_values;
  vector<double> relative_measurement_errors;

  ifstream data_file(filename);

  double Q2_settings[3] = {};

  for (int k=0; k<size(Q2_settings); k++) {
    cout << "Reading: " << filename << endl;
    string line;
    while(getline (data_file, line)) {
      long unsigned int i = 0;
      string value = "";
      while(line[i] != ' ') {
        value += line[i];
        i++;
      }
      if (stod(value) != Q2_settings[k]) {continue;}
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

  const int dim = 2;
  double res, err;

  double xl[2] = {0, 0};
  double xu[2] = {100, 1};

    struct parameters params = {1, 1};

    const gsl_rng_type *T;
    gsl_rng *rng;

    gsl_monte_function L_G = {&L_g, dim, &params};
    gsl_monte_function T_G = {&T_g, dim, &params};

    gsl_rng_env_setup();
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

      while (true) {
        status = gsl_monte_vegas_integrate(&L_G, xl, xu, dim, integration_calls, rng, L_s, &res, &err);
        if (status != 0) {throw "gsl error";}
        if (gsl_monte_vegas_chisq(L_s) < 5) {
          break;
        }
      }
      
      double sigma_L = res;

      gsl_monte_vegas_free(L_s);

      gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);

      status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
      if (status != 0) {throw "gsl error";}

      while (true) {
        status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
        if (status != 0) {throw "gsl error";}
        if (gsl_monte_vegas_chisq(T_s) < 5) {
          break;
        }
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
    stream << fixed << setprecision(0) << Q2_settings[k];
    TString title = "Reduced cross section fit for Q^{2}="+stream.str()+";x;#sigma_{r} (mb)";
    comparison_graph->SetTitle(title);

    TGraphErrors* measurement_data = new TGraphErrors(size(Q2_values), measured_x, measured_sigma, measured_x_error, measured_sigma_error);
    measurement_data->SetTitle("Reduced cross section fit");
    comparison_graph->Add(measurement_data, "P");

    TGraph* model_fit = new TGraph(size(Q2_values), measured_x, model_sigma);
    comparison_graph->Add(model_fit, "*");

    comparison_graph->Draw("A");

    TLegend* legend = new TLegend(0.65, 0.7, 0.9, 0.9);
    legend->AddEntry(measurement_data,"Measurement data");
    legend->AddEntry(model_fit,"Model fit", "P");
    legend->SetTextSize(0.04);
    legend->Draw();

    gPad->SetLogx();

    TString outfile_name = "figures/fit_comparison_Q2_"+stream.str()+".pdf";
    comparison_canvas->Print(outfile_name);

    gsl_rng_free(rng);
  }
  return 0;
}
