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

#include "cross_section_file_reader.h"


int main() {

  string numerator_filename("/home/jarno/Cross_section_calculation/data/transform_diff_L_sigma_bk_p.txt");
  string denominator_filename("/home/jarno/Cross_section_calculation/data/direct_diff_L_sigma_bk_p.txt");

  vector<double> initial_numerator_Q2, initial_denominator_Q2, initial_numerator_x, initial_numerator_sigma, initial_numerator_sigma_error, initial_denominator_x, initial_denominator_sigma, initial_denominator_sigma_error;

  read_sigma_file(numerator_filename, initial_numerator_Q2, initial_numerator_x, initial_numerator_sigma, initial_numerator_sigma_error);
  read_sigma_file(denominator_filename, initial_denominator_Q2, initial_denominator_x, initial_denominator_sigma, initial_denominator_sigma_error);

  vector<double> numerator_Q2, denominator_Q2;
  vector<vector<double>> numerator_x, numerator_sigma, numerator_sigma_error, denominator_x, denominator_sigma, denominator_sigma_error;

  split_by_Q2(numerator_Q2, numerator_x, numerator_sigma, numerator_sigma_error, initial_numerator_Q2, initial_numerator_x, initial_numerator_sigma, initial_numerator_sigma_error);
  split_by_Q2(denominator_Q2, denominator_x, denominator_sigma, denominator_sigma_error, initial_denominator_Q2, initial_denominator_x, initial_denominator_sigma, initial_denominator_sigma_error);

  TMultiGraph* comparison_graphs = new TMultiGraph();
  comparison_graphs->SetTitle("Ratio between different integration methods in longitudinal case");
  for (long unsigned int i=0; i < numerator_Q2.size(); i++) {
    double ratio[numerator_x[i].size()];
    double x[numerator_x[i].size()];
    for (long unsigned int j=0; j<numerator_x[i].size(); j++) {
      ratio[j] = numerator_sigma[i][j]/denominator_sigma[i][j];
      x[j] = numerator_x[i][j];
    }
    TGraph* subgraph = new TGraph(numerator_x[i].size(), x, ratio);
    TString subgraph_name = "Q^{2}=" + to_string(numerator_Q2[i]);
    subgraph->SetTitle(subgraph_name);
    comparison_graphs->Add(subgraph, "C");
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  if (true) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  comparison_canvas->Print("figures/L_transform_direct_ratio.pdf");
  
  return 0;
}
