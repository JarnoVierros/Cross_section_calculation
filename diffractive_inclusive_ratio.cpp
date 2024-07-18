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

  ifstream denominator_data_file("archive/data/Jdipamp/J_T_inclusive_sigma_x.txt");
  string line;
  bool first_line = true;
  while (getline (denominator_data_file, line)) {
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

  for (long unsigned int i=0; i<Q2_values.size(); i++) {
  }

  denominator_data_file.clear();
  denominator_data_file.seekg(0);

  vector<vector<double>> x_values, denominator, numerator;
  vector<double> current_x, current_denominator, current_numerator;
  int Q2_index = 0;

  first_line = true;
  while (getline (denominator_data_file, line)) {
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
      denominator.push_back(current_denominator);
      current_x = {};
      current_denominator = {};
      Q2_index++;
    }
    current_x.push_back(x_value);
    current_denominator.push_back(sigma_value);
  }
  x_values.push_back(current_x);
  denominator.push_back(current_denominator);

  ifstream numerator_data_file("archive/data/Jdipamp/diffractive/diff_T_sigma_x.txt");

  Q2_index = 0;
  current_x = {};
  first_line = true;
  while (getline (numerator_data_file, line)) {
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
      numerator.push_back(current_numerator);
      current_x = {};
      current_numerator = {};
      Q2_index++;
    }
    current_x.push_back(x_value);
    current_numerator.push_back(sigma_value);
  }
  x_values.push_back(current_x);
  numerator.push_back(current_numerator);


  TMultiGraph* comparison_graphs = new TMultiGraph();
  comparison_graphs->SetTitle("Transverse exclusive inclusive ratio");
  for (long unsigned int i=0; i < Q2_values.size(); i++) {
    double ratio[x_values[i].size()];
    double x[x_values[i].size()];
    for (long unsigned int j=0; j<x_values[i].size(); j++) {
      ratio[j] = numerator[i][j]/denominator[i][j];
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
  if (true) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  comparison_canvas->Print("figures/exclusive_inclusive_ratio.pdf");
  
  return 0;
}
