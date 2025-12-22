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
#include "TText.h"
#include "TLatex.h"

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
static double e_f = sqrt(2.0/3*2.0/3+1.0/3*1.0/3+1.0/3*1.0/3);
static double m_f = 0; //GeV

const double normalization = 8/M_PI*alpha_em*N_c*e_f*e_f;

//35, 34, 33, 32, 30, 25, 20, 15, 10, 5, 4, 3, 2, 1, 0.5
// 17, 16, 15, 10, 5, 1

static double r_limit = 34.64; // 34.64
static double b_min_limit = 17.32; // 17.32

const bool print_r_limit = false;
const bool print_b_min_limit = false;
const string dipole_amp_type = "bk";
const string nucleus_type = "p";
const string filename_end = "";

const int warmup_calls = 10000;//
const int integration_calls = 300000;//
const int integration_iterations = 1;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double shifted_x(double x, double Q2) {
  //return x;
  return x*(1+(4*m_f*m_f)/Q2);
}

double dipole_amplitude(double r, double b_min, double phi, double x) {
  if (nucleus_type == "p") {
    return get_p_dipole_amplitude(p_table, r, b_min, phi, x);
  } else if (nucleus_type == "Pb") {
    return get_Pb_dipole_amplitude(Pb_table, r, b_min, phi, x);
  } else {
    throw 1;
  }
}

double L_integrand(double r, double b_min, double phi, double z, double Q2, double x) {
  return r*b_min*4*Q2*z*z*gsl_pow_2(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r))*dipole_amplitude(r, b_min, phi, shifted_x(x, Q2));
}

double T_integrand(double r, double b_min, double phi, double z, double Q2, double x) {
  return r*b_min*(m_f*m_f*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z, Q2)*r)) + epsilon2(z, Q2)*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z, Q2)*r)))*dipole_amplitude(r, b_min, phi, shifted_x(x, Q2));
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

struct thread_par_struct
{
  double Q2;
  double x;
  double &sigma;
  double &sigma_error;
  thread_par_struct(double a1, double a2, double &a3, double &a4) : Q2(a1), x(a2), sigma(a3), sigma_error(a4) {}
};

void integrate_for_L_sigma(thread_par_struct par) {

  const int dim = 4;
  double res, err;

  double xl[dim] = {0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
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
    cout << "nan found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "L, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(L_s) << endl;

  gsl_monte_vegas_free(L_s);
}

void integrate_for_T_sigma(thread_par_struct par) {


  const int dim = 4;
  double res, err;

  double xl[dim] = {0, 0, 0, 0};
  double xu[dim] = {r_limit, b_min_limit, M_PI, 1};

  struct parameters params = {1, 1};
  params.Q2 = par.Q2;
  params.x = par.x;
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
    cout << "nan found at x=" << params.x << endl;
  }
  sigma = res;
  sigma_error = err;
  cout << "T, Q²=" << params.Q2 << ", x=" << params.x << ", res: " << sigma << ", err: " << sigma_error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);
}

int main() {

  gsl_set_error_handler_off();

  vector<double> Q2_values;
  vector<double> x_values;
  vector<double> y_values;
  vector<double> measured_sigma_values;
  vector<double> relative_measurement_errors;

  string data_filename = "data/F2_data.dat";
  ifstream data_file(data_filename);

  cout << "Reading: " << data_filename << endl;
  string line;
  while(getline (data_file, line)) {
    long unsigned int i = 0;
    
    while (line[i] != ' ') {
      i++;
    }
    i++;
    
    string value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    Q2_values.push_back(stod(value));
    i++;

    bool has_x = false;
    for (long unsigned int j=0; j<size(line); j++) {
      if (line[j] == 'x') {
        has_x = true;
        break;
      }
    }

    if (has_x) {
      string part_1 = "";
      while(line[i] != ' ') {
        part_1 += line[i];
        i++;
      }
      i+=6;
      char part_2 = line[i];
      string x_string = part_1 + "e-" + part_2;
      x_values.push_back(stod(x_string));

      i+=2;
    } else {
      value = "";
      while(line[i] != ' ') {
        value += line[i];
        i++;
      }
      x_values.push_back(stod(value));
      i++;
    }


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

    for (int j=0; j<8; j++) {
      while(line[i] != ' ') {
        i++;
      }
      i++;
    }

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    relative_measurement_errors.push_back(stod(value));
  }
  /*
  for (long unsigned int i=0; i<Q2_values.size(); i++) {
    cout << "Q2=" << Q2_values[i] << ", x=" << x_values[i] << ", y=" << y_values[i] << ", sigma=" << measured_sigma_values[i] << ", error=" << relative_measurement_errors[i] << endl;
  }
  */
  string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
  if (nucleus_type == "p") {
    load_p_dipole_amplitudes(p_table, filename);
  } else if (nucleus_type == "Pb") {
    load_Pb_dipole_amplitudes(Pb_table, filename);
  } else {
    throw 1;
  }


  vector<double> predicted_sigma_values;
  vector<double> predicted_errors;

  double predicted_L_sigmas[Q2_values.size()], predicted_T_sigmas[Q2_values.size()], predicted_L_errors[Q2_values.size()], predicted_T_errors[Q2_values.size()];
  thread L_threads[Q2_values.size()], T_threads[Q2_values.size()];

  for (long unsigned int i=0; i<Q2_values.size(); i++) {
    thread_par_struct L_par(Q2_values[i], x_values[i], predicted_L_sigmas[i], predicted_L_errors[i]);
    L_threads[i] = thread(integrate_for_L_sigma, L_par);

    thread_par_struct T_par(Q2_values[i], x_values[i], predicted_T_sigmas[i], predicted_T_errors[i]);
    T_threads[i] = thread(integrate_for_T_sigma, T_par);
  }

  for (long unsigned int i=0; i<Q2_values.size(); i++) {
    L_threads[i].join();

    T_threads[i].join();
  }

  e_f = 2.0/3;
  m_f = 1.27;

  double c_predicted_L_sigmas[Q2_values.size()], c_predicted_T_sigmas[Q2_values.size()], c_predicted_L_errors[Q2_values.size()], c_predicted_T_errors[Q2_values.size()];
  thread c_L_threads[Q2_values.size()], c_T_threads[Q2_values.size()];

  for (long unsigned int i=0; i<Q2_values.size(); i++) {
    thread_par_struct L_par(Q2_values[i], x_values[i], c_predicted_L_sigmas[i], c_predicted_L_errors[i]);
    c_L_threads[i] = thread(integrate_for_L_sigma, L_par);

    thread_par_struct T_par(Q2_values[i], x_values[i], c_predicted_T_sigmas[i], c_predicted_T_errors[i]);
    c_T_threads[i] = thread(integrate_for_T_sigma, T_par);
  }

  for (long unsigned int i=0; i<Q2_values.size(); i++) {
    c_L_threads[i].join();

    c_T_threads[i].join();
  }
  

  double predicted_sigma_r[Q2_values.size()];

  for (long unsigned int i=0; i<Q2_values.size(); i++) {
    double FL = Q2_values[i]/(4*M_PI*M_PI*alpha_em)*(predicted_L_sigmas[i] + c_predicted_L_sigmas[i]);
    double FT = Q2_values[i]/(4*M_PI*M_PI*alpha_em)*(predicted_T_sigmas[i] + c_predicted_T_sigmas[i]);
    double F2 = FL + FT;
    predicted_sigma_r[i] = F2 - y_values[i]*y_values[i]/(1+gsl_pow_2(1-y_values[i]))*FL;
  }

  vector<TMultiGraph*> comparison_graphs;
  vector<double> unique_Q2_values;
  vector<double> x, x_error, measurement, measurement_error, prediction;

  double x_limits[2] = {1e-7, 1};
  double y_limits[2] = {0, 1.5};

  double previous_Q2 = Q2_values[0];
  unique_Q2_values.push_back(Q2_values[0]);
  for (long unsigned int i=0; i<Q2_values.size(); i++) {
    if (Q2_values[i] != previous_Q2) {

      double* x_arr = &x[0];
      double* x_error_arr = &x_error[0];
      double* measurement_arr = &measurement[0];
      double* measurement_error_arr = &measurement_error[0];
      double* prediction_arr = &prediction[0];

      TGraphErrors* measurement_graph = new TGraphErrors(x.size(), x_arr, measurement_arr, x_error_arr, measurement_error_arr);

      TGraph* prediction_graph = new TGraph(x.size(), x_arr, prediction_arr);

      TMultiGraph* comparison_graph = new TMultiGraph();

      comparison_graph->Add(measurement_graph, "P");
      comparison_graph->Add(prediction_graph, "C");
      bool huh = previous_Q2 < 0.2;
      cout << "Q2=" << previous_Q2 << ", bool=" << huh << endl;
      if (previous_Q2 < 0.2) {
        comparison_graph->GetXaxis()->SetLimits(1e-7, 1e-5);
        comparison_graph->GetYaxis()->SetRangeUser(0, 0.25);
      } else {
        comparison_graph->GetXaxis()->SetLimits(x_limits[0], x_limits[1]);
        comparison_graph->GetYaxis()->SetRangeUser(y_limits[0], y_limits[1]);
      }
      
      comparison_graph->GetXaxis()->SetLabelSize(0.05);
      comparison_graph->GetYaxis()->SetLabelSize(0.05);

      //comparison_graph->GetXaxis()->SetTitle("x");
      //comparison_graph->GetYaxis()->SetTitle("#sigma_{r}");

      comparison_graphs.push_back(comparison_graph);

      x.clear();
      x_error.clear();
      measurement.clear();
      measurement_error.clear();
      prediction.clear();

      previous_Q2 = Q2_values[i];
      unique_Q2_values.push_back(Q2_values[i]);

    }
    x.push_back(x_values[i]);
    x_error.push_back(0);
    measurement.push_back(measured_sigma_values[i]);
    measurement_error.push_back(relative_measurement_errors[i]/100*measured_sigma_values[i]);
    prediction.push_back(predicted_sigma_r[i]);

  }

  double* x_arr = &x[0];
  double* x_error_arr = &x_error[0];
  double* measurement_arr = &measurement[0];
  double* measurement_error_arr = &measurement_error[0];
  double* prediction_arr = &prediction[0];

  TGraphErrors* measurement_graph = new TGraphErrors(x.size(), x_arr, measurement_arr, x_error_arr, measurement_error_arr);

  TGraph* prediction_graph = new TGraph(x.size(), x_arr, prediction_arr);

  TMultiGraph* comparison_graph = new TMultiGraph();

  comparison_graph->Add(measurement_graph, "P");
  comparison_graph->Add(prediction_graph, "C");

  comparison_graph->GetXaxis()->SetLimits(x_limits[0], x_limits[1]);
  comparison_graph->GetYaxis()->SetRangeUser(y_limits[0], y_limits[1]);
  
  comparison_graph->GetXaxis()->SetLabelSize(0.05);
  comparison_graph->GetYaxis()->SetLabelSize(0.05);

  //comparison_graph->GetXaxis()->SetTitle("x");
  //comparison_graph->GetYaxis()->SetTitle("#sigma_{r}");

  comparison_graphs.push_back(comparison_graph);


  int figure_width = 5;
  int figure_height;
  if (unique_Q2_values.size()-(unique_Q2_values.size()/5)*5 == 0) {
    figure_height = unique_Q2_values.size()/5;
  } else {
    figure_height = unique_Q2_values.size()/5 + 1;
  }
  
  TCanvas* multicanvas = new TCanvas("multicanvas", "multipads", figure_width*10000, figure_height*10000);
  multicanvas->Divide(figure_width, figure_height, 0, 0);
  cout << "size=" << comparison_graphs.size() << endl;
  for (long unsigned int i=0; i<comparison_graphs.size(); i++) {

    multicanvas->cd(i+1);

    comparison_graphs[i]->Draw("A");

    gPad->SetLogx();

    stringstream Q2_stream;
    Q2_stream << setprecision(2) << unique_Q2_values[i];
    TString Q2_string = "Q^{2}=" + Q2_stream.str();
    TLatex* Q2_text;
    if (unique_Q2_values[i] < 0.2) {
      Q2_text = new TLatex(2e-7, 0.2, Q2_string);
    } else {
      Q2_text = new TLatex(1e-6, 1.4, Q2_string);
    }
    Q2_text->Draw("Same");

    if (unique_Q2_values[i] < 0.2) {
      TLatex* axis_tick = new TLatex(8e-7, 0.02, "10^{-6}");
      axis_tick->Draw("Same");
    }

  }

  multicanvas->cd(0);

  TLatex* y_axis_label = new TLatex(0.005, 0.97, "#sigma_{r}");
  y_axis_label->SetTextSize(0.015);
  y_axis_label->Draw("Same");

  TLatex* x_axis_label = new TLatex(0.985, 0.005, "x");
  x_axis_label->SetTextSize(0.02);
  x_axis_label->Draw("Same");

  TString figure_filename = "figures/F2_data_comparison.pdf";
  multicanvas->Print(figure_filename);

}



