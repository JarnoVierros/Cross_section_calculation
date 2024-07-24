
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

#include "cross_section_file_reader.h"

const double alpha_em = 1.0/137;

void read_data_file(string filename, vector<double> &Q2_values, vector<double> &beta_values, vector<double> &x_values, vector<double> &x_pom_F2_values, vector<double> &delta_stat_values, vector<double> &delta_sys_values) {
  ifstream data_file(filename);

  cout << "Reading: " << filename << endl;
  string line;
  while(getline (data_file, line)) {

    long unsigned int i = 0;
    string value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    Q2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    beta_values.push_back(stod(value));
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
    x_pom_F2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    delta_stat_values.push_back(stod(value));
    i++;

    value = "";
    while(i < line.size()) {
      value += line[i];
      i++;
    }
    delta_sys_values.push_back(stod(value));
    i++;

  }
  cout << "Finished reading file" << endl;
}

int main() {

  double Q2_selection = 4.5;
  double beta_selection = 0.04;

  vector<double> L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit;
  string L_prediction_filenames[] = {"data/differential_diffractive_L_20mil_0-4.txt"};

  for (long unsigned int i=0; i<size(L_prediction_filenames); i++) {
    read_differential_sigma_file(L_prediction_filenames[i], L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit);
  }

  vector<double> T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit;
  string T_prediction_filenames[] = {"data/differential_diffractive_L_20mil_0-4.txt"};

  for (long unsigned int i=0; i<size(T_prediction_filenames); i++) {
    read_differential_sigma_file(T_prediction_filenames[i], T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit);
  }


  vector<double> measurement_Q2, measurement_beta, measurement_x, measurement_xpomF2, measurement_delta_stat, measurement_delta_sys;
  string measurement_filename = "data/differential_HERA_data.dat";

  read_data_file(measurement_filename, measurement_Q2, measurement_beta, measurement_x, measurement_xpomF2, measurement_delta_stat, measurement_delta_sys);

  vector<double> x_selection, chosen_measurement_xpomF2, chosen_delta;

  for (int i=0; i<measurement_Q2.size(); i++) {
    if (measurement_Q2[i] != Q2_selection) {
      continue;
    }
    if (measurement_beta[i] != beta_selection) {
      continue;
    }
    x_selection.push_back(measurement_x[i]);
    chosen_measurement_xpomF2.push_back(measurement_xpomF2[i]);
    chosen_delta.push_back(measurement_delta_stat[i]+measurement_delta_sys[i]);
  }

  vector<double> chosen_prediction_xpomF2, chosen_prediction_error;

  for (int i=0; i<x_selection.size(); i++) {
    double L_sigma, L_error, T_sigma, T_error;
    double L_found = false, T_found = false;
    for (int j=0; j<L_prediction_Q2.size(); j++) {
      if (L_prediction_Q2[i] != Q2_selection) {
        continue;
      }
      if (L_prediction_beta[i] != beta_selection) {
        continue;
      }
      if (L_prediction_x[i] != x_selection[j]) {
        continue;
      }
      L_sigma = L_prediction_sigma[i];
      L_error = L_prediction_sigma_error[i];
      L_found = true;
      break;
    }
    for (int j=0; j<T_prediction_Q2.size(); j++) {
      if (T_prediction_Q2[i] != Q2_selection) {
        continue;
      }
      if (T_prediction_beta[i] != beta_selection) {
        continue;
      }
      if (T_prediction_x[i] != x_selection[j]) {
        continue;
      }
      T_sigma = T_prediction_sigma[i];
      T_error = T_prediction_sigma_error[i];
      T_found = true;
      break;
    }
    if (!L_found) {
      cout << "Warning: L prediction not found" << endl;
      L_sigma = 0;
    }
    if (!T_found) {
      cout << "Warning: T prediction not found" << endl;
      T_sigma = 0;
    }
    chosen_prediction_xpomF2[i] = Q2_selection*Q2_selection/(pow(2*M_PI, 2)*alpha_em*beta_selection)*(L_sigma + T_sigma);
    chosen_prediction_error[i] = Q2_selection*Q2_selection/(pow(2*M_PI, 2)*beta_selection*alpha_em)*sqrt(L_error*L_error + T_error*T_error);
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);

  double zeroes[x_selection.size()];
  for (int i=0; i<x_selection.size(); i++) {
    zeroes[i] = 0;
  }
  double* x_selection_arr = &x_selection[0];
  double* chosen_measurement_xpomF2_arr = &chosen_measurement_xpomF2[0];
  double* chosen_delta_arr = &chosen_delta[0];

  TGraphErrors* measurement_data = new TGraphErrors(x_selection.size(), x_selection_arr, chosen_measurement_xpomF2_arr, zeroes, chosen_delta_arr);
  measurement_data->SetTitle("x_pom*F_2^D(3) data prediction comparison;x;x_pom*F_2^D(3)");

  double* chosen_prediction_xpomF2_arr = &chosen_prediction_xpomF2[0];
  double* chosen_prediction_error_arr = &chosen_prediction_error[0];

  TGraphErrors* prediction = new TGraphErrors(x_selection.size(), x_selection_arr, chosen_prediction_xpomF2_arr, zeroes, chosen_prediction_error_arr);
  prediction->SetLineColor(2);

  float location[4];
  location[0] = 0.65;
  location[1] = 0.7;
  location[2] = 0.9;
  location[3] = 0.9;

  TLegend* legend = new TLegend(location[0], location[1], location[2], location[3]);
  legend->AddEntry(measurement_data,"Measurement data");
  legend->AddEntry(prediction,"Prediction");
  legend->SetTextSize(0.04);
  legend->Draw();

  comparison_canvas->Print("figures/F2D_data_comparison.pdf");

  return 0;
}
