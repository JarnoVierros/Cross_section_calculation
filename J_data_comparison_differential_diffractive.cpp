
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

struct plot {
  TMultiGraph* comparison_graph;
  TGraphErrors* measurement_data;
  TGraphErrors* prediction;
  TGraph* FL_prediction;
  TGraph* FT_prediction;
};

int main() {

  vector<plot> plots;

  vector<double> L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit;
  string L_prediction_filenames[] = {"data/differential_diffractive_L_20mil_0-4.txt", "data/differential_diffractive_L_20mil_5-19.txt", "data/differential_diffractive_L_20mil_20-39.txt"};

  for (long unsigned int i=0; i<size(L_prediction_filenames); i++) {
    read_differential_sigma_file(L_prediction_filenames[i], L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit);
  }

  /*
  for (int i=0; i<L_prediction_Q2.size(); i++) {
    cout << L_prediction_x[i] << endl;
  }
  */

  vector<double> T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit;
  string T_prediction_filenames[] = {"data/differential_diffractive_T_20mil_0-4.txt", "data/differential_diffractive_T_20mil_5-19.txt", "data/differential_diffractive_T_20mil_20-39.txt"};

  for (long unsigned int i=0; i<size(T_prediction_filenames); i++) {
    read_differential_sigma_file(T_prediction_filenames[i], T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit);
  }


  vector<double> measurement_Q2, measurement_beta, measurement_x, measurement_xpomF2, measurement_delta_stat, measurement_delta_sys;
  string measurement_filename = "data/differential_HERA_data.dat";

  read_data_file(measurement_filename, measurement_Q2, measurement_beta, measurement_x, measurement_xpomF2, measurement_delta_stat, measurement_delta_sys);


  double Q2_selections[] = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 7.5, 7.5, 7.5, 7.5};
  double beta_selections[] = {0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 0.04, 0.1, 0.2, 0.4};



  for (int k=0; k<size(Q2_selections); k++) {

    vector<double> x_selection, chosen_measurement_xpomF2, chosen_delta;

    for (int i=0; i<measurement_Q2.size(); i++) {
      if (measurement_Q2[i] != Q2_selections[k]) {
        continue;
      }
      if (measurement_beta[i] != beta_selections[k]) {
        continue;
      }
      x_selection.push_back(measurement_x[i]);
      chosen_measurement_xpomF2.push_back(measurement_xpomF2[i]);
      chosen_delta.push_back(measurement_delta_stat[i]+measurement_delta_sys[i]);
    }


    vector<double> chosen_prediction_xpomF2, chosen_prediction_xpomFL, chosen_prediction_xpomFT, chosen_prediction_error;

    for (int i=0; i<x_selection.size(); i++) {
      double L_sigma, L_error, T_sigma, T_error;
      double L_found = false, T_found = false;
      for (int j=0; j<L_prediction_Q2.size(); j++) {

        if (L_prediction_Q2[j] != Q2_selections[k]) {
          continue;
        }
        if (L_prediction_beta[j] != beta_selections[k]) {
          continue;
        }
        if (L_prediction_x[j] != x_selection[i]) {
          continue;
        }
        L_sigma = L_prediction_sigma[j];
        L_error = L_prediction_sigma_error[j];
        L_found = true;
        break;
      }
      for (int j=0; j<T_prediction_Q2.size(); j++) {
        if (T_prediction_Q2[j] != Q2_selections[k]) {
          continue;
        }
        if (T_prediction_beta[j] != beta_selections[k]) {
          continue;
        }
        if (T_prediction_x[j] != x_selection[i]) {
          continue;
        }
        T_sigma = T_prediction_sigma[j];
        T_error = T_prediction_sigma_error[j];
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
      chosen_prediction_xpomF2.push_back(Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*(L_sigma + T_sigma));
      chosen_prediction_xpomFL.push_back(Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*L_sigma);
      chosen_prediction_xpomFT.push_back(Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*T_sigma);
      chosen_prediction_error.push_back(Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*beta_selections[k]*alpha_em)*sqrt(L_error*L_error + T_error*T_error));
    }
/*
struct plot {
  TString title;
  TGraphErrors* measurement_data;
  TGraphErrors* prediction;
  TGraph* FL_prediction;
  TGraph* FT_prediction;
};
*/

    double zeroes[x_selection.size()];
    for (int i=0; i<x_selection.size(); i++) {
      zeroes[i] = 0;
    }
    double* x_selection_arr = &x_selection[0];
    double* chosen_measurement_xpomF2_arr = &chosen_measurement_xpomF2[0];
    double* chosen_delta_arr = &chosen_delta[0];

    TMultiGraph* comparison_graph = new TMultiGraph();
    stringstream Q2_stream;
    Q2_stream << fixed << setprecision(1) << Q2_selections[k];
    stringstream beta_stream;
    beta_stream << fixed << setprecision(2) << beta_selections[k];
    TString title = "x_{pom}F_{2}^{D(3)} measurement prediction comparison at Q^{2}="+Q2_stream.str()+", #beta="+beta_stream.str()+";x;x_{pom}F_{2}^{D(3)}";
    comparison_graph->SetTitle(title);

    TGraphErrors* measurement_data = new TGraphErrors(x_selection.size(), x_selection_arr, chosen_measurement_xpomF2_arr, zeroes, chosen_delta_arr);
    measurement_data->SetTitle("Measurement data");
    comparison_graph->Add(measurement_data, "P");

    double* chosen_prediction_xpomF2_arr = &chosen_prediction_xpomF2[0];
    double* chosen_prediction_error_arr = &chosen_prediction_error[0];

    TGraphErrors* prediction = new TGraphErrors(x_selection.size(), x_selection_arr, chosen_prediction_xpomF2_arr, zeroes, chosen_prediction_error_arr);
    prediction->SetTitle("Prediction");
    prediction->SetMarkerStyle(8);
    prediction->SetMarkerColor(2);
    prediction->SetLineColor(2);
    comparison_graph->Add(prediction, "C");

    double* chosen_prediction_xpomFL_arr = &chosen_prediction_xpomFL[0];

    TGraph* FL_prediction = new TGraph(x_selection.size(), x_selection_arr, chosen_prediction_xpomFL_arr);
    FL_prediction->SetLineColor(3);
    comparison_graph->Add(FL_prediction, "C");

    double* chosen_prediction_xpomFT_arr = &chosen_prediction_xpomFT[0];

    TGraph* FT_prediction = new TGraph(x_selection.size(), x_selection_arr, chosen_prediction_xpomFT_arr);
    FT_prediction->SetLineColor(4);
    comparison_graph->Add(FT_prediction, "C");

    TString title = "x_{pom}F_{2}^{D(3)} measurement prediction comparison at Q^{2}="+Q2_stream.str()+", #beta="+beta_stream.str()+";x;x_{pom}F_{2}^{D(3)}";
    plot new_plot = {comparison_graph, measurement_data, prediction, FL_prediction, FT_prediction};
    plots.push_back(new_plot);

  }

  TCanvas* multicanvas = new TCanvas("multicanvas", "multipads", 900, 700);
  multicanvas->Divide(6, 8);

  for (int i=0; i<plots.size(); i++) {
    multicanvas->cd(i);
    plots[i].comparison_graph->Draw("A");

    float location[4];
    location[0] = 0.5;
    location[1] = 0.7;
    location[2] = 0.75;
    location[3] = 0.9;

    TLegend* legend = new TLegend(location[0], location[1], location[2], location[3]);
    legend->AddEntry(plots[i].measurement_data,"Measurement data");
    legend->AddEntry(plots[i].prediction,"Prediction");
    legend->AddEntry(plots[i].FL_prediction,"F_{L}");
    legend->AddEntry(plots[i].FT_prediction,"F_{T}");
    legend->SetTextSize(0.04);
    legend->Draw();
  }

  TString figure_filename = "figures/F2D_data_comparison_Q2.pdf";
  multicanvas->Print(figure_filename);

  return 0;
}
