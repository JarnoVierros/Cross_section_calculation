#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
using namespace std;


int main() {

  vector<int> Q2_values;

  ifstream inclusive_data_file("data/inclusive_T_sigma_x.txt");
  string line;
  bool first_line = true;
  while (getline (inclusive_data_file, line)) {
    if (first_line) {
      first_line = false;
      continue;
    }
    string Q2_string = "";
    int i = 0;
    while (line[i] != ';') {
      Q2_string += line[i];
      i++;
    }
    double Q2_value = stod(Q2_string);
    bool found = false;
    for (long unsigned int i=0; i<Q2_values.size(); i++) {
      if (Q2_value == Q2_values[i]) {
        found = true;
        break;
      }
    }
    if (!found) {
      Q2_values.push_back(Q2_value);
    }
  }

  for (int i=0; i<Q2_values.size(); i++) {
  }

  inclusive_data_file.clear();
  inclusive_data_file.seekg(0);

  vector<vector<double>> x_values, inclusive_sigma, diff_sigma;
  vector<double> current_x, current_inclusive_sigma, current_diff_sigma;
  int Q2_index = 0;

  first_line = true;
  while (getline (inclusive_data_file, line)) {
    if (first_line) {
      first_line = false;
      continue;
    }
    string Q2_string = "";
    int i = 0;
    while (line[i] != ';') {
      Q2_string += line[i];
      i++;
    }
    i++;
    double Q2_value = stod(Q2_string);

    string x_value_string = "";
    while (line[i] != ';') {
      x_value_string += line[i];
      i++;
    }
    i++;
    double x_value = stod(x_value_string);
    string sigma_value_string = "";
    while (line[i] != ';') {
      sigma_value_string += line[i];
      i++;
    }
    i++;
    double sigma_value = stod(sigma_value_string);
    if (Q2_values[Q2_index] != Q2_value) {
      x_values.push_back(current_x);
      inclusive_sigma.push_back(current_inclusive_sigma);
      current_x = {};
      current_inclusive_sigma = {};
      Q2_index++;
    }
    current_x.push_back(x_value);
    current_inclusive_sigma.push_back(sigma_value);
  }
  x_values.push_back(current_x);
  inclusive_sigma.push_back(current_inclusive_sigma);

  ifstream diffractive_data_file("data/diff_T_sigma_x.txt");

  Q2_index = 0;
  current_x = {};
  first_line = true;
  while (getline (diffractive_data_file, line)) {
    if (first_line) {
      first_line = false;
      continue;
    }
    string Q2_string = "";
    int i = 0;
    while (line[i] != ';') {
      Q2_string += line[i];
      i++;
    }
    i++;
    double Q2_value = stod(Q2_string);

    string x_value_string = "";
    while (line[i] != ';') {
      x_value_string += line[i];
      i++;
    }
    i++;
    double x_value = stod(x_value_string);

    string sigma_value_string = "";
    while (line[i] != ';') {
      sigma_value_string += line[i];
      i++;
    }
    i++;
    double sigma_value = stod(sigma_value_string);
    if (Q2_values[Q2_index] != Q2_value) {
      x_values.push_back(current_x);
      diff_sigma.push_back(current_diff_sigma);
      current_x = {};
      current_inclusive_sigma = {};
      Q2_index++;
    }
    current_x.push_back(x_value);
    current_diff_sigma.push_back(sigma_value);
  }
  x_values.push_back(current_x);
  diff_sigma.push_back(current_diff_sigma);


  TMultiGraph* comparison_graphs = new TMultiGraph();
  comparison_graphs->SetTitle("Ratio between diffractive and inclusive transverse cross sections");
  for (long unsigned int i=0; i < Q2_values.size(); i++) {
    double ratio[x_values[i].size()];
    double x[x_values[i].size()];
    for (int j=0; j<x_values[i].size(); j++) {
      ratio[j] = diff_sigma[i][j]/inclusive_sigma[i][j];
      x[j] = x_values[i][j];
    }
    TGraph* subgraph = new TGraph(x_values[i].size(), x, ratio);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[i]);
    subgraph->SetTitle(subgraph_name);
    comparison_graphs->Add(subgraph);
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graphs->Draw("A PMC PLC");

  gPad->SetLogx();

  comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);

  comparison_canvas->Print("figures/diffractive_inclusive_T_ratio.pdf");
  
  return 0;
}